#!/bin/bash

# Parse input arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bam_file_list) BAM_FILE_LIST="$2"; shift ;;
        --gtf_file) GTF_FILE="$2"; shift ;;
        --splice_table) SPLICE_TABLE="$2"; shift ;;
        --output_dir) OUTPUT_DIR="$2"; shift ;;
        --n_threads) N_THREADS="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Check if required arguments are provided
if [[ -z "$BAM_FILE_LIST" || -z "$GTF_FILE" || -z "$SPLICE_TABLE" || -z "$OUTPUT_DIR" || -z "$N_THREADS" ]]; then
    echo "Usage: $0 --bam_file_list <path_to_bam_file_list> --gtf_file <path_to_gtf_file> --splice_table <path_to_splice_table> --output_dir <output_directory> --n_threads <number_of_threads>"
    exit 1
fi

# Create the output directory
mkdir -p "$OUTPUT_DIR"

# Run LeafCutter
echo "Running LeafCutter..."
python src/run_leafcutter.py --bam_file_list "$BAM_FILE_LIST" --output_dir "$OUTPUT_DIR" --n_threads "$N_THREADS"

# Process LeafCutter PSI
echo "Processing LeafCutter PSI..."
python src/process_leafcutter_psi.py --input_file "$OUTPUT_DIR/leafcutter/leafcutter_psi.txt.gz" --output_file "$OUTPUT_DIR/leafcutter/leafcutter_per_site_psi.csv"

# Add LeafCutter PSI values to splice table
echo "Adding LeafCutter PSI values to splice table..."
python src/leafcutter_to_splice_table.py --leafcutter_file "$OUTPUT_DIR/leafcutter/leafcutter_per_site_psi.csv" --splice_table_file "$SPLICE_TABLE" --output_file "$OUTPUT_DIR/leafcutter/leafcutter_final_splice_table.txt" --n_threads "$N_THREADS"

# Run rMATS
echo "Running rMATS..."
python src/run_rmats.py --bam_file_list "$BAM_FILE_LIST" --output_dir "$OUTPUT_DIR" --gtf_file "$GTF_FILE" --n_threads "$N_THREADS"

# Process rMATS PSI
echo "Processing rMATS PSI..."
python src/process_rmats_psi.py --input_dir "$OUTPUT_DIR/rmats/results" --outfile "$OUTPUT_DIR/rmats/rmats_per_site_psi.csv" --n_threads "$N_THREADS"

# Add rMATS PSI values to splice table
echo "Adding rMATS PSI values to splice table..."
python src/rmats_to_splice_table.py --rmats_file "$OUTPUT_DIR/rmats/rmats_per_site_psi.csv" --splice_table_file "$SPLICE_TABLE" --output_file "$OUTPUT_DIR/rmats/rmats_final_splice_table.txt"

# Run SpliSER
echo "Running SpliSER..."
python src/run_spliser.py --bam_file_list "$BAM_FILE_LIST" --output_dir "$OUTPUT_DIR" --gtf_file "$GTF_FILE" --n_threads "$N_THREADS"

# Add SpliSER values to splice table
echo "Adding SpliSER values to splice table..."
python src/spliser_to_splice_table.py --spliser_output "$OUTPUT_DIR/spliser/spliser_outputAll.DiffSpliSER.tsv" --splice_table "$SPLICE_TABLE" --output_file "$OUTPUT_DIR/spliser/spliser_final_splice_table.txt" --n_threads "$N_THREADS"

echo "All steps completed successfully!"


# usage:
# ./src/run_all.sh --bam_file_list test_data/bam_files.tsv --gtf_file GRCh38/gencode.v29.primary_assembly.annotation.gtf.gz --splice_table output/GRCh38_v29_splice_table.txt --output_dir output --n_threads 20
