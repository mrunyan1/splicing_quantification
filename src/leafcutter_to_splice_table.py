import pandas as pd
import argparse
from concurrent.futures import ThreadPoolExecutor

def process_row(idx, row, df):
    current_chromosome = row[2]  # Assuming chromosome is in the second column
    tx_start = row[4]  # Assuming tx_start is in the 4th column
    tx_end = row[5]  # Assuming tx_end is in the 5th column

    # Find rows in df where chromosome matches and position is within tx_start and tx_end
    matching_rows = df[(df['chromosome'] == current_chromosome) & (df['position'] >= tx_start) & (df['position'] <= tx_end)]

    # Extract positions and PSI values, join them as comma-separated strings with a trailing comma
    positions = ','.join(map(str, matching_rows['position'].tolist())) + ','
    psi_values = ','.join(map(str, matching_rows['psi'].tolist())) + ','

    return idx, positions, psi_values

def process_splice_table(leafcutter_file, splice_table_file, output_file, n_threads):
    # Load the input files
    df = pd.read_csv(leafcutter_file)
    splice_table = pd.read_csv(splice_table_file, sep="\t", header=None)

    # Add new columns to splice_table
    splice_table['positions'] = ''  # Initialize as empty string
    splice_table['PSI'] = ''  # Initialize as empty string

    splice_table[4] = pd.to_numeric(splice_table[4], errors='coerce')  # tx_start
    splice_table[5] = pd.to_numeric(splice_table[5], errors='coerce')  # tx_end

    # Use ThreadPoolExecutor for multithreading
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = {executor.submit(process_row, idx, row, df): idx for idx, row in splice_table.iterrows()}

        for future in futures:
            idx, positions, psi_values = future.result()
            splice_table.at[idx, 'positions'] = positions
            splice_table.at[idx, 'PSI'] = psi_values

    # Remove rows where 'positions' or 'PSI' columns contain only a comma
    splice_table = splice_table[~splice_table['positions'].isin([',']) & ~splice_table['PSI'].isin([','])]

    # Save the modified splice_table to the output file
    splice_table.to_csv(output_file, sep="\t", index=False, header=False)

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Process splice table with leafcutter PSI data using multithreading.")
    parser.add_argument('--leafcutter_file', type=str, required=True, help="Path to the input leafcutter PSI CSV file.")
    parser.add_argument('--splice_table_file', type=str, required=True, help="Path to the input splice table file.")
    parser.add_argument('--output_file', type=str, required=True, help="Path to the output file.")
    parser.add_argument('--n_threads', type=int, default=4, help="Number of threads to use for parallel processing (default: 4)")

    args = parser.parse_args()

    # Run the processing function
    process_splice_table(args.leafcutter_file, args.splice_table_file, args.output_file, args.n_threads)

if __name__ == "__main__":
    main()
    
# usage:
# python leafcutter_to_splice_table.py --leafcutter_file /path/to/leafcutter_per_site_psi.csv --splice_table_file /path/to/splice_table.txt --output_file /path/to/output_file.txt --n_threads 8

