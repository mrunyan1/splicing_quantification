import os
import subprocess
import logging
import pandas as pd
import argparse
import gzip
import shutil

def decompress_gtf(gtf_file, output_dir):
    """
    Decompress the gzipped GTF file if necessary and return the path to the decompressed file.
    """
    if gtf_file.endswith('.gz'):
        logging.info(f"GTF file {gtf_file} is gzipped. Decompressing...")
        decompressed_gtf_file = os.path.join(output_dir, os.path.basename(gtf_file).replace('.gz', ''))
        with gzip.open(gtf_file, 'rb') as f_in, open(decompressed_gtf_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        logging.info(f"Decompressed GTF file to {decompressed_gtf_file}")
        return decompressed_gtf_file
    else:
        return gtf_file

def collect_bam_files_from_tsv(tsv_file):
    """
    Collect BAM files from the TSV file and return a dictionary
    mapping sample names to lists of BAM files.
    """
    bam_pairs = {}
    
    try:
        # Read the TSV file into a DataFrame
        bam_df = pd.read_csv(tsv_file, delimiter='\t', header=None)
        bam_file_list_dir = os.path.dirname(tsv_file)  # Get the directory of the TSV file

        # Each row represents a sample; the columns contain paths to BAM files for the replicates
        for index, row in bam_df.iterrows():
            # Use the base name of the first BAM file (without the extension) as the sample name
            first_bam_base = os.path.splitext(os.path.basename(row[0]))[0]
            bam_files = []

            for file in row:
                if isinstance(file, str) and file.endswith(".bam"):
                    bam_path = os.path.join(bam_file_list_dir, file)
                    if os.path.exists(bam_path):
                        bam_files.append(bam_path)
                    else:
                        logging.warning(f"BAM file not found: {bam_path}")

            if len(bam_files) >= 2:  # Ensure there are at least two BAM files
                bam_pairs[first_bam_base] = bam_files
            else:
                logging.warning(f"Expected at least 2 BAM files, but found {len(bam_files)} for sample {first_bam_base}")

    except Exception as e:
        logging.error(f"Error reading BAM files from TSV: {e}")
    
    logging.info(f"Total samples with paired BAM files: {len(bam_pairs)}")
    
    return bam_pairs


def run_rmats(bam_files, output_dir, gtf_file, sample_name, n_threads):
    """
    Function to run rMATS using a list of BAM files
    """
    output_rmats = os.path.join(output_dir, "rmats", "results", sample_name)
    
    if os.path.isdir(output_rmats):
        summary_file = os.path.join(output_rmats, "summary.txt")
        if os.path.isfile(summary_file):
            logging.info(f"rMATS output already exists for sample {sample_name}. Skipping.")
            return
        else:
            logging.info(f"Output directory exists but summary.txt is missing for sample {sample_name}. Continuing with processing.")

    os.makedirs(output_rmats, exist_ok=True)
    tmp_dir = os.path.join(output_rmats, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    bam_list_file = os.path.join(output_rmats, "bam_list.txt")
    with open(bam_list_file, 'w') as f:
        f.write(','.join(bam_files) + '\n')

    command = [
        "rmats.py",
        "--b1", bam_list_file,
        "--gtf", gtf_file,
        "--od", output_rmats,
        "--tmp", tmp_dir,
        "--statoff",
        "--readLength", "100",
        "--variable-read-length",
        "--nthread", str(n_threads),
        "-t", "paired",
        "--allow-clipping",
        "--libType", "fr-firststrand",
        "--novelSS"
    ]

    logging.info(f"Running rMATS with BAM files: {bam_files}")
    subprocess.run(command, check=True)
    logging.info(f"rMATS finished, results are in {output_rmats}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run rMATS on BAM files listed in a TSV file.")
    parser.add_argument("--bam_file_list", required=True, help="Path to the TSV file listing BAM files.")
    parser.add_argument("--output_dir", required=True, help="Path to the output directory.")
    parser.add_argument("--gtf_file", required=True, help="Path to the GTF file (gzipped or uncompressed).")
    parser.add_argument("--n_threads", type=int, default=54, help="Number of threads to use for rMATS (default: 54).")
    args = parser.parse_args()
    
    bam_file_list = args.bam_file_list
    output_dir = args.output_dir
    gtf_file = args.gtf_file
    n_threads = args.n_threads

    # Set up logging
    log_dir = os.path.join(output_dir, "rmats")
    os.makedirs(log_dir, exist_ok=True)
    log_file_path = os.path.join(log_dir, "process_log.txt")
    logging.basicConfig(filename=log_file_path, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    logging.info("Processing started.")
    
    # Decompress the GTF file if it's gzipped
    gtf_file = decompress_gtf(gtf_file, output_dir)

    # Collect BAM files from the TSV file
    bam_pairs = collect_bam_files_from_tsv(bam_file_list)
    
    processed_samples = 0
    total_samples = len(bam_pairs)
    
    for sample_name, bam_files in bam_pairs.items():
        run_rmats(bam_files, output_dir, gtf_file, sample_name, n_threads)
        processed_samples += 1
        logging.info(f"Processed {processed_samples} out of {total_samples} samples.")
    
    logging.info("Processing finished.")

if __name__ == "__main__":
    main()

# usage:
# python run_rmats.py --bam_file_list test_data/bam_files.tsv --output_dir output --gtf_file GRCh38/gencode.v29.primary_assembly.annotation.gtf.gz --n_threads 20
