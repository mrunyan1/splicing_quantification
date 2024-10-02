import pandas as pd
import argparse

def process_leafcutter_data(input_file, output_file):
    df = pd.read_csv(input_file, compression='gzip', sep=None, engine='python')

    df = df.reset_index()
    df = df.rename(columns={'index': 'intron_info'})

    df['intron_info'] = df['intron_info'].astype(str)

    # split the 'intron_info' column into separate columns
    df[['chromosome', 'intron_start', 'intron_end', 'cluster_id_strand']] = df['intron_info'].str.split(':', expand=True)

    df['strand'] = df['cluster_id_strand'].str[-1]  # Last character is the strand (+ or -)
    df['cluster_id'] = df['cluster_id_strand'].str[:-1]  # Everything except the last character
    df['cluster_id'] = df['cluster_id'].str.rstrip('_')  # Remove trailing underscore from cluster_id

    df['intron_start'] = df['intron_start'].astype(int)
    df['intron_end'] = df['intron_end'].astype(int)


    df = df.drop(columns=['cluster_id_strand'])

    # exon end/start sites by subtracting/adding 1 from intron start/end
    df['exon_end'] = df['intron_start'] - 1
    df['exon_start'] = df['intron_end'] + 1

    # calculate the average of all sample PSI values for each row, ignoring NA values
    sample_columns = [col for col in df.columns if col not in ['intron_info', 'chromosome', 'intron_start', 'intron_end', 'exon_start', 'exon_end', 'cluster_id', 'strand']]
    df['average_psi'] = df[sample_columns].mean(axis=1, skipna=True)

    # reorder
    columns_order = ['intron_info', 'chromosome', 'intron_start', 'intron_end', 'exon_start', 'exon_end', 'cluster_id', 'strand', 'average_psi'] + \
                    [col for col in df.columns if col not in ['intron_info', 'chromosome', 'intron_start', 'intron_end', 'exon_start', 'exon_end', 'cluster_id', 'strand', 'average_psi']]

    df = df[columns_order]

    # dict to accumulate average psi values for each unique position within clusters
    position_accumulator = {}

    # iter over each row in the df
    for _, row in df.iterrows():
        cluster_id = row['cluster_id']
        intron_info = row['intron_info']
        chromosome = row['chromosome']
        strand = row['strand']
        
        exon_start = row['exon_start']
        exon_end = row['exon_end']
        
        average_psi = row['average_psi']
        
        # add exon_start position to dict
        if (cluster_id, exon_start, 'exon_start') not in position_accumulator:
            position_accumulator[(cluster_id, exon_start, 'exon_start')] = {
                'psi': 0,
                'intron_info': intron_info,
                'chromosome': chromosome,
                'strand': strand
            }
        position_accumulator[(cluster_id, exon_start, 'exon_start')]['psi'] += average_psi
        
        # add exon_end position to dict
        if (cluster_id, exon_end, 'exon_end') not in position_accumulator:
            position_accumulator[(cluster_id, exon_end, 'exon_end')] = {
                'psi': 0,
                'intron_info': intron_info,
                'chromosome': chromosome,
                'strand': strand
            }
        position_accumulator[(cluster_id, exon_end, 'exon_end')]['psi'] += average_psi

    # create dataframe from dict
    summed_position_df = pd.DataFrame([
        {
            'cluster_id': key[0], 
            'position': key[1], 
            'type': key[2], 
            'psi': value['psi'], 
            'intron_info': value['intron_info'],
            'chromosome': value['chromosome'],
            'strand': value['strand']
        } 
        for key, value in position_accumulator.items()
    ])

    # Save the new dataframe to output file
    summed_position_df.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser(description="Process LeafCutter PSI data.")
    parser.add_argument('--input_file', type=str, required=True, help='Path to the input gzip-compressed PSI file.')
    parser.add_argument('--output_file', type=str, required=True, help='Path to the output CSV file.')

    args = parser.parse_args()

    process_leafcutter_data(args.input_file, args.output_file)

if __name__ == "__main__":
    main()

# usage:
# python process_leafcutter_psi.py --input_file /path/to/input/leafcutter_psi.txt.gz --output_file /path/to/output/leafcutter_per_site_psi.csv
