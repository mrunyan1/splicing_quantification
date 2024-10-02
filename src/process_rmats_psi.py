"""
This script processes rMATS output files from multiple directories to extract splice site information and calculate average Percent Spliced In (PSI) values for splice sites across the 5 event types for one results directory. It identifies donor and acceptor splice sites specific to each event type, organizes them into DataFrames, and combines the results into a single matrix where splice sites are rows and directories are columns. 
"""
import pandas as pd
import glob
import os
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed


def extract_splice_sites(event_type, file_path):
    """
    Extracts splice site information from rMATS output files for a specified splicing event type.

    This function reads a rMATS output file corresponding to a specific splicing event type 
    (e.g., SE, MXE, A3SS, A5SS, RI) and extracts the splice sites for each event in the file.
    
    Inclusion and Skip Forms by Event Type:
    - SE (Skipped Exon):
      - Inclusion Form: Includes the exon between 'exonStart_0base' and 'exonEnd'.
      - Skip Form: Skips the exon, so no donor or acceptor sites from the exon are included.

    - MXE (Mutually Exclusive Exons):
      - Inclusion Form:
        - If strand is '+', includes the first exon between '1stExonStart_0base' and '1stExonEnd'.
        - If strand is '-', includes the second exon between '2ndExonStart_0base' and '2ndExonEnd'.
      - Skip Form:
        - If strand is '+', skips the first exon and includes the second exon.
        - If strand is '-', skips the second exon and includes the first exon.

    - A5SS (Alternative 5' Splice Site):
      - Inclusion Form: Includes the long exon defined by 'longExonStart_0base' and 'longExonEnd'.
      - Skip Form: 
        - If strand is '-', uses 'shortES' as the skip site.
        - If strand is '+', uses 'shortEE' as the skip site.
      - PSI Handling: 
        - If strand is '-', the `longExonEnd` PSI is always 1 because it is included in both forms.
        - If strand is '+', the `longExonStart_0base` PSI is always 1 because it is included in both forms.

    - A3SS (Alternative 3' Splice Site):
      - Inclusion Form: Includes the long exon defined by 'longExonStart_0base' and 'longExonEnd'.
      - Skip Form: 
        - If strand is '-', uses 'shortEE' as the skip site.
        - If strand is '+', uses 'shortES' as the skip site.
      - PSI Handling:
        - If strand is '-', the `longExonStart_0base` PSI is always 1 because it is included in both forms.
        - If strand is '+', the `longExonEnd` PSI is always 1 because it is included in both forms.

    - RI (Retained Intron):
      - Inclusion Form: Retains the intron, so the donor and acceptor sites are 'riExonStart_0base' and 'riExonEnd'.
      - Skip Form: The intron is not retained, so uses 'upstreamEE' and 'downstreamES' as donor and acceptor sites.

    Parameters:
    event_type (str): The type of splicing event (e.g., 'SE', 'MXE', 'A3SS', 'A5SS', 'RI').
    file_path (str): The file path to the rMATS output file containing the splicing event data.

    Returns:
    pd.DataFrame: A DataFrame containing the extracted splice site information, including columns 
    for splice site coordinates, PSI values, event type, and other relevant details.
    """
    df = pd.read_csv(file_path, sep='\t')

    def safe_psi_calculation(psi_column):
        def calculate(psi):
            values = psi.split(',')
            valid_values = [float(v) for v in values if v != 'NA']
            return sum(valid_values) / len(valid_values) if valid_values else None
        return df[psi_column].apply(calculate)

    df['psi_1'] = safe_psi_calculation('IncLevel1')
    df['skip_psi_1'] = df['psi_1'].apply(lambda x: 1 - x if x is not None else None)

    result_df = pd.DataFrame()

    if event_type == 'SE':
        donor_sites = df['exonStart_0base']
        acceptor_sites = df['exonEnd']

        donor_df = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': donor_sites,
            'site_type': 'donor',
            'event_type': event_type,
            'form_type': 'inclusion',
            'psi_1': df['psi_1'],
            'FDR': df['FDR']
        })

        acceptor_df = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': acceptor_sites,
            'site_type': 'acceptor',
            'event_type': event_type,
            'form_type': 'inclusion',
            'psi_1': df['psi_1'],
            'FDR': df['FDR']
        })

        result_df = pd.concat([donor_df, acceptor_df], ignore_index=True)

    elif event_type == 'MXE':
        donor_sites_incl = df.apply(lambda x: x['1stExonStart_0base'] if x['strand'] == '+' else x['2ndExonStart_0base'], axis=1)
        acceptor_sites_incl = df.apply(lambda x: x['1stExonEnd'] if x['strand'] == '+' else x['2ndExonEnd'], axis=1)
        donor_sites_skip = df.apply(lambda x: x['2ndExonStart_0base'] if x['strand'] == '+' else x['1stExonStart_0base'], axis=1)
        acceptor_sites_skip = df.apply(lambda x: x['2ndExonEnd'] if x['strand'] == '+' else x['1stExonEnd'], axis=1)

        donor_df_incl = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': donor_sites_incl,
            'site_type': 'donor',
            'event_type': event_type,
            'form_type': 'inclusion',
            'psi_1': df['psi_1'],
            'FDR': df['FDR']
        })

        acceptor_df_incl = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': acceptor_sites_incl,
            'site_type': 'acceptor',
            'event_type': event_type,
            'form_type': 'inclusion',
            'psi_1': df['psi_1'],
            'FDR': df['FDR']
        })

        donor_df_skip = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': donor_sites_skip,
            'site_type': 'donor',
            'event_type': event_type,
            'form_type': 'skip',
            'psi_1': df['skip_psi_1'],
            'FDR': df['FDR']
        })

        acceptor_df_skip = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': acceptor_sites_skip,
            'site_type': 'acceptor',
            'event_type': event_type,
            'form_type': 'skip',
            'psi_1': df['skip_psi_1'],
            'FDR': df['FDR']
        })

        result_df = pd.concat([donor_df_incl, acceptor_df_incl, donor_df_skip, acceptor_df_skip], ignore_index=True)

    if event_type == 'A5SS':
        donor_sites = df['longExonStart_0base']
        acceptor_sites_incl = df['longExonEnd']
        acceptor_sites_skip = df.apply(lambda x: x['shortES'] if x['strand'] == '-' else x['shortEE'], axis=1)

        donor_df = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': donor_sites,
            'site_type': 'donor',
            'event_type': event_type,
            'form_type': 'inclusion',
            'psi_1': df.apply(lambda x: 1 if x['strand'] == '+' else x['psi_1'], axis=1),  # PSI always 1 for longExonStart if strand is +
            'FDR': df['FDR']
        })

        acceptor_df_incl = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': acceptor_sites_incl,
            'site_type': 'acceptor',
            'event_type': event_type,
            'form_type': 'inclusion',
            'psi_1': df.apply(lambda x: 1 if x['strand'] == '-' else x['psi_1'], axis=1),  # PSI always 1 for longExonEnd if strand is -
            'FDR': df['FDR']
        })

        acceptor_df_skip = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': acceptor_sites_skip,
            'site_type': 'acceptor',
            'event_type': event_type,
            'form_type': 'skip',
            'psi_1': df['skip_psi_1'],
            'FDR': df['FDR']
        })

        result_df = pd.concat([donor_df, acceptor_df_incl, acceptor_df_skip], ignore_index=True)

    elif event_type == 'A3SS':
        donor_sites = df['longExonStart_0base']
        acceptor_sites_incl = df['longExonEnd']
        acceptor_sites_skip = df.apply(lambda x: x['shortEE'] if x['strand'] == '-' else x['shortES'], axis=1)

        donor_df = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': donor_sites,
            'site_type': 'donor',
            'event_type': event_type,
            'form_type': 'inclusion',
            'psi_1': df.apply(lambda x: 1 if x['strand'] == '-' else x['psi_1'], axis=1),  # PSI always 1 for longExonStart if strand is -
            'FDR': df['FDR']
        })

        acceptor_df_incl = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': acceptor_sites_incl,
            'site_type': 'acceptor',
            'event_type': event_type,
            'form_type': 'inclusion',
            'psi_1': df.apply(lambda x: 1 if x['strand'] == '+' else x['psi_1'], axis=1),  # PSI always 1 for longExonEnd if strand is +
            'FDR': df['FDR']
        })

        acceptor_df_skip = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': acceptor_sites_skip,
            'site_type': 'acceptor',
            'event_type': event_type,
            'form_type': 'skip',
            'psi_1': df['skip_psi_1'],
            'FDR': df['FDR']
        })

        result_df = pd.concat([donor_df, acceptor_df_incl, acceptor_df_skip], ignore_index=True)

    elif event_type == 'RI':
        donor_sites_incl = df['riExonStart_0base']
        acceptor_sites_incl = df['riExonEnd']
        donor_sites_skip = df['upstreamEE']
        acceptor_sites_skip = df['downstreamES']

        donor_df_incl = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': donor_sites_incl,
            'site_type': 'donor',
            'event_type': event_type,
            'form_type': 'inclusion',
            'psi_1': df['psi_1'],
            'FDR': df['FDR']
        })

        acceptor_df_incl = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': acceptor_sites_incl,
            'site_type': 'acceptor',
            'event_type': event_type,
            'form_type': 'inclusion',
            'psi_1': df['psi_1'],
            'FDR': df['FDR']
        })

        donor_df_skip = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': donor_sites_skip,
            'site_type': 'donor',
            'event_type': event_type,
            'form_type': 'skip',
            'psi_1': df['skip_psi_1'],
            'FDR': df['FDR']
        })

        acceptor_df_skip = pd.DataFrame({
            'ID': df['ID'],
            'geneSymbol': df['geneSymbol'],
            'chr': df['chr'],
            'strand': df['strand'],
            'splice_site': acceptor_sites_skip,
            'site_type': 'acceptor',
            'event_type': event_type,
            'form_type': 'skip',
            'psi_1': df['skip_psi_1'],
            'FDR': df['FDR']
        })

        result_df = pd.concat([donor_df_incl, acceptor_df_incl, donor_df_skip, acceptor_df_skip], ignore_index=True)

    return result_df




def process_directory(rmats_dir, event_files):
    """
    Processes a single rMATS results directory to extract and average PSI values.
    """
    all_splice_sites = []

    # check for required event files in the directory
    valid_directory = any(os.path.exists(os.path.join(rmats_dir, file_name)) for file_name in event_files.values())
    if not valid_directory:
        print(f"Skipping {rmats_dir} as it does not contain any required event files.")
        return None

    # iter over each event type and corresponding file name
    for event_type, file_name in event_files.items():
        file_path = os.path.join(rmats_dir, file_name)
        if os.path.exists(file_path):
            splice_sites_df = extract_splice_sites(event_type, file_path)
            all_splice_sites.append(splice_sites_df)
        else:
            print(f"File not found for {event_type} in {rmats_dir}: {file_path}")

    if all_splice_sites:
        all_splice_sites_df = pd.concat(all_splice_sites, ignore_index=True)
        
        # group by splice_site, chr, and geneSymbol to calculate average PSI
        average_psi_df = all_splice_sites_df.groupby(['splice_site', 'chr', 'geneSymbol']).agg(
            average_psi=('psi_1', 'mean')
        ).reset_index()

        # rename the column for the directory name
        dir_name = os.path.basename(os.path.normpath(rmats_dir))
        average_psi_df = average_psi_df.rename(columns={'average_psi': dir_name})

        # set the splice site, chromosome, and gene name as index
        average_psi_df.set_index(['splice_site', 'chr', 'geneSymbol'], inplace=True)

        return average_psi_df
    else:
        print(f"No valid splice site data found in {rmats_dir}.")
        return None


def main():
    parser = argparse.ArgumentParser(description="Process rMATS data to calculate average PSI across samples.")
    parser.add_argument('--input_dir', required=True, help='Path to the parent directory containing rMATS result directories.')
    parser.add_argument('--outfile', required=True, help='Path to the output CSV file for average PSI across samples.')
    parser.add_argument('--n_threads', type=int, default=20, help='Number of threads to use for processing (default: 20)')

    args = parser.parse_args()

    rmats_parent_dir = args.input_dir
    output_file = args.outfile
    n_threads = args.n_threads

    rmats_dirs = [os.path.join(rmats_parent_dir, d) for d in os.listdir(rmats_parent_dir) 
                  if os.path.isdir(os.path.join(rmats_parent_dir, d)) and not d.startswith('.')]

    event_files = {
        'SE': 'SE.MATS.JC.txt',
        'MXE': 'MXE.MATS.JC.txt',
        'A3SS': 'A3SS.MATS.JC.txt',
        'A5SS': 'A5SS.MATS.JC.txt',
        'RI': 'RI.MATS.JC.txt'
    }

    all_average_psi = {}

    # ThreadPoolExecutor for multithreading with specified number of threads
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        # submit a processing task for each rMATS directory
        futures = {executor.submit(process_directory, rmats_dir, event_files): rmats_dir for rmats_dir in rmats_dirs}

        for future in as_completed(futures):
            dir_name = os.path.basename(os.path.normpath(futures[future]))
            try:
                result = future.result()
                if result is not None:
                    all_average_psi[dir_name] = result
            except Exception as e:
                print(f"Error processing {dir_name}: {e}")

    # combine all average PSI data into a df if there are results
    if all_average_psi:
        combined_average_psi_df = pd.concat(all_average_psi.values(), axis=1)

        # calculate the average PSI across all samples for each splice site
        average_across_samples = combined_average_psi_df.mean(axis=1, skipna=True)
        
        # combine with the index (splice_site, chr, geneSymbol) from the combined df
        average_across_samples_df = combined_average_psi_df.reset_index()[['splice_site', 'chr', 'geneSymbol']].copy()
        average_across_samples_df['average_psi_across_samples'] = average_across_samples.values
        
        average_across_samples_df = average_across_samples_df.dropna(subset=['average_psi_across_samples'])
        
        average_across_samples_df.to_csv(output_file, index=False)
        print(f"Average PSI across samples saved to {output_file}")


if __name__ == "__main__":
    main()
    
# usage:
# python process_rmats_psi.py --input_dir /path/to/rmats_directory --outfile /path/to/output_file.csv --threads 20
