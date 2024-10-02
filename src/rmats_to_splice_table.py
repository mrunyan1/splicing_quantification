import pandas as pd
import argparse

def process_rmats_splice_table(rmats_file, splice_table_file, output_file):

    df = pd.read_csv(rmats_file)
    splice_table = pd.read_csv(splice_table_file, sep="\t", header=None)

    # add new columns to splice_table
    splice_table['positions'] = ''  # Initialize as an empty string
    splice_table['PSI'] = ''  # Initialize as an empty string

    # Ensure tx_start and tx_end are numeric
    splice_table[4] = pd.to_numeric(splice_table[4], errors='coerce')  # tx_start
    splice_table[5] = pd.to_numeric(splice_table[5], errors='coerce')  # tx_end

    # Iterate over each row in the splice_table
    for idx, row in splice_table.iterrows():
        current_chromosome = row[2]  # Assuming chromosome is in the second column
        tx_start = row[4]  # Assuming tx_start is in the 4th column
        tx_end = row[5]  # Assuming tx_end is in the 5th column

        # Find rows in df where chromosome matches and position is within tx_start and tx_end
        matching_rows = df[(df['chr'] == current_chromosome) & (df['splice_site'] >= tx_start) & (df['splice_site'] <= tx_end)]

        # Extract positions and PSI values, join them as comma-separated strings with a trailing comma
        positions = ','.join(map(str, matching_rows['splice_site'].tolist())) + ','
        psi_values = ','.join(map(str, matching_rows['average_psi_across_samples'].tolist())) + ','

        # Update the splice_table row with the concatenated positions and PSI values
        splice_table.at[idx, 'positions'] = positions
        splice_table.at[idx, 'PSI'] = psi_values

    # Remove rows where 'positions' or 'PSI' columns contain only a comma
    splice_table = splice_table[~splice_table['positions'].isin([',']) & ~splice_table['PSI'].isin([','])]

    # Save the modified splice_table to the output file
    splice_table.to_csv(output_file, sep="\t", index=False, header=False)

def main():

    parser = argparse.ArgumentParser(description="Process splice table with rMATS PSI data.")
    parser.add_argument('--rmats_file', type=str, required=True, help="Path to the input rMATS PSI CSV file.")
    parser.add_argument('--splice_table_file', type=str, required=True, help="Path to the input splice table file.")
    parser.add_argument('--output_file', type=str, required=True, help="Path to the output file.")

    args = parser.parse_args()


    process_rmats_splice_table(args.rmats_file, args.splice_table_file, args.output_file)

if __name__ == "__main__":
    main()
    
# usage:
# python rmats_to_splice_table.py --rmats_file /path/to/rmats_psi.csv --splice_table_file /path/to/splice_table.txt --output_file /path/to/output_file.txt

