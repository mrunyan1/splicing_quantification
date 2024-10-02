import sys
import argparse
import pandas as pd
import numpy as np
import gzip
import re
import os

def extract_attribute(attr_string, attr_name):
    pattern = f'{attr_name} "([^"]+)"'
    match = re.search(pattern, attr_string)
    return match.group(1) if match else None

def main():
    parser = argparse.ArgumentParser(description="Process GTF file and check paralog status using pre-downloaded paralog data.")
    parser.add_argument('--gtf_file', type=str, required=True, help='Path to the input GTF file (gzipped).')
    parser.add_argument('--output_file', type=str, required=True, help='Path to the output file.')
    args = parser.parse_args()

    gtf_file = args.gtf_file
    paralog_file = "test_data/paralogs_GRCh38.txt"
    output_file = args.output_file

    print("Reading GTF file...")

    with gzip.open(gtf_file, 'rt') as f:
        gtf_data = pd.read_csv(f, sep='\t', comment='#', header=None, 
                               names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"])

    gtf_data_parsed = gtf_data.copy()
    for attr in ['gene_id', 'transcript_id', 'gene_name', 'gene_type', 'transcript_type', 'level', 'tag', 'havana_gene', 'havana_transcript']:
        gtf_data_parsed[attr] = gtf_data_parsed['attributes'].apply(lambda x: extract_attribute(x, attr))
    gtf_data_parsed = gtf_data_parsed.drop('attributes', axis=1)

    # protein-coding genes
    protein_coding_genes = gtf_data_parsed[gtf_data_parsed['gene_type'] == 'protein_coding']

    # identify primary transcripts using appris tags if available
    primary_transcripts = protein_coding_genes[protein_coding_genes['feature'] == 'transcript'].copy()
    primary_transcripts['appris_priority'] = primary_transcripts['tag'].apply(lambda x: 
        1 if 'appris_principal_1' in str(x) else
        2 if 'appris_principal_2' in str(x) else
        3 if 'appris_principal_3' in str(x) else
        4 if 'appris_principal_4' in str(x) else
        5 if 'appris_principal_5' in str(x) else 6)

    primary_transcripts = primary_transcripts.sort_values(['gene_id', 'appris_priority', 'end', 'start'], 
                                                          ascending=[True, True, False, True])
    primary_transcripts = primary_transcripts.groupby('gene_id').first().reset_index()
    primary_transcripts = primary_transcripts[['gene_id', 'transcript_id']].rename(columns={'transcript_id': 'primary_transcript_id'})

    principal_transcripts_data = protein_coding_genes.merge(primary_transcripts, 
                                                            left_on='transcript_id', 
                                                            right_on='primary_transcript_id')

    exon_sites = gtf_data_parsed[gtf_data_parsed['feature'] == 'exon'].groupby('transcript_id').agg({
        'start': lambda x: ','.join(map(str, sorted(x))),
        'end': lambda x: ','.join(map(str, sorted(x)))
    }).reset_index()
    exon_sites.columns = ['transcript_id', 'exon_start_sites', 'exon_end_sites']

    # Merge exon sites and filter
    principal_transcripts_with_exons = principal_transcripts_data.merge(exon_sites, on='transcript_id')
    principal_transcripts_with_exons = principal_transcripts_with_exons[principal_transcripts_with_exons['feature'] == 'transcript']
    principal_transcripts_with_exons_filtered = principal_transcripts_with_exons[
        principal_transcripts_with_exons['exon_start_sites'].str.count(',') >= 1
    ].copy()  # Explicitly create a copy here

    # Trim the exon start and end sites to include a trailing comma
    principal_transcripts_with_exons_filtered['exon_start_sites'] = principal_transcripts_with_exons_filtered['exon_start_sites'].apply(lambda x: ','.join(x.split(',')[1:]) + ',')
    principal_transcripts_with_exons_filtered['exon_end_sites'] = principal_transcripts_with_exons_filtered['exon_end_sites'].apply(lambda x: ','.join(x.split(',')[:-1]) + ',')

    # Use the filtered DataFrame
    principal_transcripts_with_exons_final = principal_transcripts_with_exons_filtered[['seqname', 'start', 'end', 'strand', 'transcript_id', 'gene_name', 'gene_id_x', 'exon_start_sites', 'exon_end_sites']]

    print("Reading pre-downloaded paralog data...")
    paralog_data = pd.read_csv(paralog_file, sep='\t')

    # Filter paralog data to rows where 'Human paralogue gene stable ID' is not empty
    paralog_data_filtered = paralog_data[paralog_data['Human paralogue gene stable ID'].notna()]

    # Paralog status based on 'gene_id_x'
    final_df = principal_transcripts_with_exons_final.copy()
    final_df['paralog_status'] = final_df['gene_id_x'].isin(paralog_data_filtered['Gene stable ID version']).astype(int)

    # Order the final DataFrame
    final_df_ordered = final_df[['transcript_id', 'paralog_status', 'seqname', 'strand', 'start', 'end', 'exon_end_sites', 'exon_start_sites']]
    
    # Create output directory if it does not exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Save to the output file
    final_df_ordered.to_csv(output_file, sep='\t', index=False, header=False)

    print(f"File saved as {output_file}")

if __name__ == "__main__":
    main()
