import os
import subprocess
import logging
import argparse
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
import pysam
import shutil
import glob

write_lock = Lock()

def extract_junctions(bam_file, junction_files_dir, junction_list_file):
    logging.info(f"Extracting junctions from {bam_file}...")
    
    # Generate junction file from BAM using regtools
    junc_file = os.path.join(junction_files_dir, os.path.basename(bam_file).replace(".bam", ".junc"))

    # Check if the junction file already exists
    if os.path.exists(junc_file):
        logging.info(f"Junction file {junc_file} already exists. Skipping extraction.")
        return
    
    try:
        pysam.index(bam_file)
    except Exception as e:
        logging.error(f"Failed to index BAM file {bam_file}: {e}")
        return
    
    # Then, extract junctions
    command_junc = [
        "regtools", "junctions", "extract", 
        "-a", "8",            # anchor length
        "-m", "20",           # Minimum intron size
        "-M", "1000000",       # Maximum intron size
        "-s", "FR",           # Strandness mode for forward-stranded data
        bam_file,
        "-o", junc_file       # Output junction file
    ]
    subprocess.run(command_junc, check=True)

    # Write the junction file to the junction list
    with write_lock:
        with open(junction_list_file, 'a') as junc_list:
            junc_list.write(junc_file + '\n')
    
    logging.info(f"Junction extraction finished, results are in {junc_file}")


def cluster_introns(junction_list_file, leafcutter_dir):
    logging.info("Clustering introns using LeafCutter...")

    # Create an output directory for the clustering files
    clustering_output = os.path.join(leafcutter_dir, "leafcutter_clustering")

    # Create a new directory for per-sample clusters
    per_sample_clusters_dir = os.path.join(leafcutter_dir, "per_sample_clusters")
    os.makedirs(per_sample_clusters_dir, exist_ok=True)

    # Run the intron clustering using LeafCutter's regtools clustering script
    command_cluster = [
        "python", "./leafcutter/clustering/leafcutter_cluster_regtools.py",
        "-j", junction_list_file,
        "-m", "1",
        "-o", clustering_output,  
        "-l", "1000000"
    ]

    
    subprocess.run(command_cluster, check=True)

    # move individual clustering files to the per_sample_clusters directory
    for file in glob.glob("*.junc.leafcutter_clustering.sorted.gz"):
        shutil.move(file, os.path.join(per_sample_clusters_dir, os.path.basename(file)))

    logging.info(f"Intron clustering finished. Overall clustering results are in {clustering_output}, per-sample files in {per_sample_clusters_dir}")

def quantify_psi(psi_output_file, leafcutter_dir):
    logging.info("Running PSI quantification...")
    
    # Assume counts_file is the output from the clustering step
    counts_file = os.path.join(leafcutter_dir, "leafcutter_clustering_perind_numers.counts.gz")

    # Run the PSI quantification command
    command_psi = [
        "Rscript", "./leafcutter/scripts/leafcutter_quantify_psi.R", 
        counts_file,
        "--output_file", psi_output_file
    ]

    subprocess.run(command_psi, check=True)
    logging.info(f"PSI quantification finished, results are in {psi_output_file}")

def main():
    parser = argparse.ArgumentParser(description="LeafCutter pipeline with PSI quantification")
    parser.add_argument('--run_psi_only', action='store_true', help="Run only the PSI quantification step")
    parser.add_argument('--bam_file_list', type=str, required=True, help="TSV file containing BAM file paths (multiple columns for replicates)")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to store the output")
    parser.add_argument('--n_threads', type=int, default=50, help="Number of threads for parallel processing (default: 50)")
    args = parser.parse_args()

    # define directories for output
    output_dir = args.output_dir
    leafcutter_dir = os.path.join(output_dir, "leafcutter")
    junction_files_dir = os.path.join(leafcutter_dir, "junction_files")
    junction_list_file = os.path.join(leafcutter_dir, "junction_files.txt")
    psi_output_file = os.path.join(leafcutter_dir, "leafcutter_psi.txt.gz")


    os.makedirs(junction_files_dir, exist_ok=True)  


    log_file = os.path.join(leafcutter_dir, "leafcutter_pipeline.log")
    logging.basicConfig(filename=log_file, filemode='a', format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

    # Initialize the junction list file
    if not os.path.exists(junction_list_file):
        with open(junction_list_file, 'w') as f:
            pass

    logging.info("LeafCutter pipeline started")
    
    if args.run_psi_only:
        # run only PSI quantification
        quantify_psi(psi_output_file, leafcutter_dir)
    else:
        # extract junctions, cluster introns, and quantify PSI
        bam_files = []

        # get the directory of the bam_file_list
        bam_file_list_dir = os.path.dirname(args.bam_file_list)

        # read the TSV file containing bam file paths
        with open(args.bam_file_list, 'r') as f:
            for line in f:
                # split the line into BAM files (handles arbitrary number of columns)
                bam_replicates = line.strip().split("\t")
                #  full paths to each BAM file and add to the list
                for bam in bam_replicates:
                    bam_files.append(os.path.join(bam_file_list_dir, bam))

        # ThreadPoolExecutor to parallelize the extraction of junctions
        max_workers = args.n_threads 
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            executor.map(lambda bam: extract_junctions(bam, junction_files_dir, junction_list_file), bam_files)

        # clustering after all junctions are extracted
        cluster_introns(junction_list_file, leafcutter_dir)
        
        # PSI quantification
        quantify_psi(psi_output_file, leafcutter_dir)
    
    logging.info("LeafCutter pipeline finished")


if __name__ == "__main__":
    main()

# usage:
# python run_leafcutter.py --bam_file_list bam_files.tsv --output_dir /path/to/output --n_threads 20
