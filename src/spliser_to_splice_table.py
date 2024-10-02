import pandas as pd
import argparse
from concurrent.futures import ThreadPoolExecutor

def process_row(row, mean_sse_df_clean):
    chrom = row.iloc[2]  # Use chromosome information
    tx_start = row.iloc[4]
    tx_end = row.iloc[5]

    # Filter the mean_sse_df_clean for the same chromosome and within transcription start/end
    filtered_sse = mean_sse_df_clean[
        (mean_sse_df_clean["Region"] == chrom) & 
        (mean_sse_df_clean["Site"] >= tx_start) & 
        (mean_sse_df_clean["Site"] <= tx_end)
    ]

    # Get the splice sites and mean SSE values
    splice_sites = filtered_sse["Site"].tolist()
    SSE_values = filtered_sse["mean_SSE"].tolist()

    # Create a dictionary with the updated values for this row
    return {
        "splice_sites": ",".join(map(str, splice_sites)) + ",",
        "SSE_values": ",".join(map(str, SSE_values)) + ","
    }

def main():
    parser = argparse.ArgumentParser(description="Process splicing data and calculate mean SSE values.")
    parser.add_argument('--spliser_output', required=True, help='Path to the input spliser output file.')
    parser.add_argument('--splice_table', required=True, help='Path to the input splice table file.')
    parser.add_argument('--output_file', required=True, help='Path to the output file for the updated splice table.')
    parser.add_argument('--n_threads', type=int, default=4, help='Number of threads for parallel processing.')

    args = parser.parse_args()

    spliser_output = pd.read_csv(args.spliser_output, sep='\t')

    # Calculate mean SSE, alpha, and beta values
    sse_columns = [col for col in spliser_output.columns if col.endswith("SSE")]
    alpha_columns = [col for col in spliser_output.columns if col.endswith("alpha")]
    beta_columns = [col for col in spliser_output.columns if col.endswith("beta")]

    spliser_output['mean_SSE'] = spliser_output[sse_columns].mean(axis=1, skipna=True)
    spliser_output['mean_alpha'] = spliser_output[alpha_columns].mean(axis=1, skipna=True)
    spliser_output['mean_beta'] = spliser_output[beta_columns].mean(axis=1, skipna=True)

    selected_columns = ['Region', 'Site', 'Strand', 'Gene', 'mean_SSE', 'mean_alpha', 'mean_beta']
    mean_sse_df = spliser_output[selected_columns]

    mean_sse_df_clean = mean_sse_df.dropna(subset=['mean_SSE'])

    splice_table = pd.read_csv(args.splice_table, sep='\t')

    # New columns for storing splice sites and SSE values
    splice_table["splice_sites"] = ""
    splice_table["SSE_values"] = ""

    # Use ThreadPoolExecutor to process rows in parallel
    with ThreadPoolExecutor(max_workers=args.n_threads) as executor:
        futures = {executor.submit(process_row, row, mean_sse_df_clean): i for i, row in splice_table.iterrows()}

        for future in futures:
            i = futures[future]
            result = future.result()

            # Update the splice table with results
            splice_table.at[i, "splice_sites"] = result["splice_sites"]
            splice_table.at[i, "SSE_values"] = result["SSE_values"]

    # Remove rows where both 'splice_sites' and 'SSE_values' are empty
    filtered_splice_table = splice_table[(splice_table['splice_sites'] != ',') | (splice_table['SSE_values'] != ',')]

    # Save the updated splice table
    filtered_splice_table.to_csv(args.output_file, sep='\t', index=False, header=None)
    print(f"Splice table saved to {args.output_file}")

if __name__ == "__main__":
    main()
