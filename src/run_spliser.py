import os
import subprocess
import logging
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
import csv

write_lock = Lock()  # lock for writing to the samples file

def create_spliser_junctions(bam_file, junction_files_dir):
    sample_name = os.path.basename(bam_file).replace(".bam", "")
    junc_file = os.path.join(junction_files_dir, f"{sample_name}.spliser.junc")
    
    if os.path.exists(junc_file):
        logging.info(f"Junction file {junc_file} already exists, skipping creation.")
        return junc_file
    
    logging.info(f"Creating junction file {junc_file} from {bam_file}...")
    
    command_junc = [
        "regtools", "junctions", "extract",
        "-a", "8", "-m", "50", "-M", "500000", "-s", "FR",
        bam_file, "-o", junc_file
    ]
    subprocess.run(command_junc, check=True)
    
    logging.info(f"Junction file created at {junc_file}")
    return junc_file


def run_spliser_process(bam_file, junc_file, process_output_dir, gtf_file):
    sample_name = os.path.basename(bam_file).replace(".bam", "")
    output_tsv = os.path.join(process_output_dir, f"{sample_name}")
    
    # Skip if the TSV file with .SpliSER.tsv suffix already exists
    if os.path.exists(f"{output_tsv}.SpliSER.tsv"):
        logging.info(f"SpliSER TSV {output_tsv}.SpliSER.tsv already exists, skipping processing.")
        return f"{output_tsv}.SpliSER.tsv"
    
    logging.info(f"Running SpliSER process on {bam_file}...")
    command_process = [
        "python", "./SpliSER/SpliSER_v0_1_8.py", "process",
        "-B", bam_file, "-b", junc_file, "-o", output_tsv,
        "-A", gtf_file, "--isStranded", "-s", "fr", "-m", "500000"
    ]
    subprocess.run(command_process, check=True)
    
    logging.info(f"SpliSER process completed, output saved to {output_tsv}")
    return output_tsv

def process_bam_file(bam_file, junction_files_dir, process_output_dir, gtf_file, samples_file):
    sample_name = os.path.basename(bam_file).replace(".bam", "")
    
    junc_file = create_spliser_junctions(bam_file, junction_files_dir)
    output_tsv = run_spliser_process(bam_file, junc_file, process_output_dir, gtf_file)
    
    if junc_file and output_tsv:
        # Safely append to the samples file
        with write_lock:
            logging.info(f"Appending {sample_name} to the samples file.")
            with open(samples_file, 'a') as f:
                f.write(f"{sample_name}\t{output_tsv}.SpliSER.tsv\t{bam_file}\n")
            logging.info(f"Successfully appended {sample_name} to the samples file.")

def run_spliser_combine(samples_file, output_combined, n_threads):
    logging.info(f"Running SpliSER combine...")
    os.makedirs(output_combined, exist_ok=True)

    command_combine = [
        "python", "./SpliSER/SpliSER_v0_1_8.py", "combine",
        "-S", samples_file, "-o", output_combined, "--isStranded", "-s", "fr", "--n_jobs", str(n_threads)
    ]
    subprocess.run(command_combine, check=True)
    logging.info(f"SpliSER combine finished, results are in {output_combined}/combined.tsv")

def run_spliser_output(samples_file, combined_file, output_diffspliser):
    logging.info(f"Running SpliSER output...")
    os.makedirs(output_diffspliser, exist_ok=True)

    command_output = [
        "python", "./SpliSER/SpliSER_v0_1_8.py", "output",
        "-S", samples_file, "-C", combined_file, "-t", "DiffSpliSER", "-o", output_diffspliser
    ]
    subprocess.run(command_output, check=True)
    logging.info(f"SpliSER output finished, results are in {output_diffspliser}")

def get_bam_files_from_tsv(tsv_file):
    """
    Get full paths for BAM files listed in the TSV file.
    """
    bam_files = []
    bam_file_list_dir = os.path.dirname(tsv_file)  # Get the directory of the TSV file
    
    with open(tsv_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) > 0:
                # Construct the full path for each BAM file in the row
                for bam in row:
                    full_bam_path = os.path.join(bam_file_list_dir, bam)
                    bam_files.append(full_bam_path)
                    
    return bam_files

def main():
    parser = argparse.ArgumentParser(description="Run SpliSER pipeline with optional steps.")
    parser.add_argument('--skip_processing', action='store_true', help="Skip processing BAM files and move directly to combining.")
    parser.add_argument('--bam_file_list', type=str, required=True, help="Path to the TSV file containing the list of BAM files.")
    parser.add_argument('--output_dir', type=str, required=True, help="Path to the output directory.")
    parser.add_argument('--gtf_file', type=str, required=True, help="Path to the GTF file.")
    parser.add_argument('--n_threads', type=int, default=12, help="Number of threads for parallel processing (default: 12)")
    args = parser.parse_args()

    spliser_dir = os.path.join(args.output_dir, "spliser")
    junction_files_dir = os.path.join(spliser_dir, "junction_files")
    samples_file = os.path.join(spliser_dir, "samples_file.tsv")
    process_output_dir = os.path.join(spliser_dir, "process")
    output_combined = os.path.join(spliser_dir, "spliser_combined")
    output_diffspliser = os.path.join(spliser_dir, "spliser_output")
    
    os.makedirs(spliser_dir, exist_ok=True)
    os.makedirs(junction_files_dir, exist_ok=True)
    os.makedirs(process_output_dir, exist_ok=True)
    os.makedirs(output_combined, exist_ok=True)
    os.makedirs(output_diffspliser, exist_ok=True)

    log_file = os.path.join(spliser_dir, "spliser_pipeline.log")
    logging.basicConfig(filename=log_file, filemode='a', format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

    logging.info("SpliSER pipeline started")
    
    if not args.skip_processing:
        bam_files_to_process = get_bam_files_from_tsv(args.bam_file_list)
        
        num_bam_files = len(bam_files_to_process)
        logging.info(f"Number of BAM files to process: {num_bam_files}")

        if not os.path.exists(samples_file):
            logging.info(f"Creating samples file at {samples_file}")
            with open(samples_file, 'w') as f:
                pass

        with ThreadPoolExecutor(max_workers=args.n_threads) as executor:
            futures = [executor.submit(process_bam_file, bam_file, junction_files_dir, process_output_dir, args.gtf_file, samples_file) for bam_file in bam_files_to_process]
            for future in as_completed(futures):
                future.result()

    run_spliser_combine(samples_file, output_combined, args.n_threads)
    combined_file = output_combined + ".combined.tsv"
    run_spliser_output(samples_file, combined_file, output_diffspliser)
    
    logging.info("SpliSER pipeline finished")

if __name__ == "__main__":
    main()
    
# usage:
# python run_spliser.py --bam_file_list ../test_data/bam_files.tsv --output_dir ../output --gtf_file ../GRCh38/gencode.v29.primary_assembly.annotation.gtf --n_threads 24