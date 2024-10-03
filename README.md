
This repository contains scripts to run the splicing quantification tools Leafcutter, rMATS and SpliSER. The output from rMATS and Leafcutter is processed to convert PSI (percent splice-in) values exons or introns to splice-site PSI values. The results from each tool are also formatted as a "splice table" which maps each splice-site (and PSI value) to its respective transcript.

## Setup
1. Clone this repository and move into the project directory:

    ```
   git clone https://github.com/mrunyan1/splicing_quantification.git
    ```

    ```
   cd splicing_quantification
    ```

3. Run the following make command to make the conda environment and download required repos:
    
    ```
   make setup
    ```

    ```
   conda activate splicing
    ```

    ```
   make install_renv
    ```

5. (workshop) From within the project directory use this command to get the test data:

   ```
   cp -r /work/talisman/mrunyan/test_data .
   ```




## Usage
### Required data
- BAM files: RNA-seq data in BAM format.
- BAM file list: A TSV file containing paths to BAM files. The TSV should be in the same location as the BAM files
and the paths should be relative to the TSV file location. Each line can contain multiple columns representing replicate BAM files.
- GTF file: A GTF file containing gene annotations. The GTF file can be downloaded using the `getReference.py` script.


### General Usage
1. Run the `getReference.py` script to download the GTF.

```
python src/getReference.py --version 29 --download gtf --output_dir GRCh38
```

- --version: Specify the version number of the GTF file to download (e.g., 29).
- --download: Choose whether to download GTF and/or reference genome (gtf, ref, or both).
- --output_dir: Directory where the downloaded files will be saved.

2. Run `make_splice_table.py` script to generate a table with principal transcript information. This table will be populated with splice-site PSI values and respective coordinates.

```
python src/make_splice_table.py --gtf_file GRCh38/gencode.v29.primary_assembly.annotation.gtf.gz --output_file output/GRCh38_v29_splice_table.txt
```

- --gtf_file: Path to the GTF file.
- --output_file: Path to the output file.

3. To run all tools/scripts that are detailed below use the `run_all.sh` script. This script will run Leafcutter, rMATS, and
SpliSER and process the output from each tool.
- Run:
```
./src/run_all.sh --bam_file_list test_data/bam_files.tsv \
                 --gtf_file GRCh38/gencode.v29.primary_assembly.annotation.gtf.gz \
                 --splice_table output/GRCh38_v29_splice_table.txt \
                 --output_dir output --n_threads 20
```

### Leafcutter
`run_leafcutter.py`

- Description: This script extracts splice junctions from BAM files using regtools, 
clusters the junctions into introns using LeafCutter, and calculates PSI
values for each intron cluster.
- Run:
  ```
  python src/run_leafcutter.py --bam_file_list test_data/bam_files.tsv --output_dir output --n_threads 20
  ```
  - `--bam_file_list`: TSV file containing BAM file paths.
  - `--output_dir`: Directory where output files will be saved.
  - `--n_threads`: Number of threads to use.
- Output:
  - Junction Files: Stored in `output_dir/leafcutter/junction_files`
  - Intron Clustering: Results saved in `output_dir/leafcutter`
  - PSI Quantification: The computed PSI values are saved in `output_dir/leafcutter/leafcutter_psi.txt.gz`

`process_leafcutter_psi.py`

- Description: Processes the PSI data from LeafCutter. Converts intron start/end sites to exon start/end sites 
and calculates PSI values for each exon start/end site.
- Run:
  ```python src/process_leafcutter_psi.py --input_file output/leafcutter/leafcutter_psi.txt.gz --output_file output/leafcutter/leafcutter_per_site_psi.csv
  ```
    - `--input_file`: Path to the LeafCutter PSI file.
    - `--output_file`: Path to the output file.
- Output: 
  - A CSV file containing the chromosome, position, PSI values, and additional information for each exon start/end site.

`leafcutter_to_splice_table.py`

- Description: For each trasncript in the splice table this scripts adds respective PSI values and genomic positions from LeafCutter quantification.
- Run:
  ```
  python src/leafcutter_to_splice_table.py \
        --leafcutter_file output/leafcutter/leafcutter_per_site_psi.csv \
        --splice_table_file output/GRCh38_v29_splice_table.txt \
        --output_file output/leafcutter/leafcutter_GRCh38_v29_splice_table.txt \
        --n_threads 20
  ```
    - `--leafcutter_file`: Path to the LeafCutter per site PSI file.
    - `--splice_table_file`: Path to the splice table file.
    - `--output_file`: Path to the output file.
    - `--n_threads`: Number of threads to use.
- Output:
  - A splice table with PSI values and positions from LeafCutter added to the table.


### rMATS
`run_rmats.py`

- Description: Runs rMATS using the BAM files provided from the bam files TSV and the GTF file.
- Run:
  ```
  python src/run_rmats.py --bam_file_list test_data/bam_files.tsv --output_dir output --gtf_file GRCh38/gencode.v29.primary_assembly.annotation.gtf.gz --n_threads 20
  ```
    - `--bam_file_list`: TSV file containing BAM file paths.
    - `--output_dir`: Directory where output files will be saved.
    - `--gtf_file`: Path to the GTF file.
    - `--n_threads`: Number of threads to use.
- Output:
  - rMATS output files are saved in `output_dir/rmats/results` see rMATS documentation for more information on the output files.


`process_rmats_psi.py`

- Description: Processes the exon PSI (inclusion level) values from rMATs to obtain site-level PSI values.
- Run:
  ```
  python src/process_rmats_psi.py --input_dir output/rmats/results --outfile output/rmats/rmats_per_site_psi.csv --n_threads 20
  ```
    - `--input_dir`: Directory containing rMATS output files.
    - `--outfile`: Path to the output file.
    - `--n_threads`: Number of threads to use.
- Output:
  - A CSV file containing the chromosome, position, PSI values, and additional information for each exon start/end site.


`rmats_to_splice_table.py`

- Description: For each trasncript in the splice table this scripts adds respective PSI values and genomic positions from rMATS quantification.
- Run:
  ```
  python src/rmats_to_splice_table.py \
      --rmats_file output/rmats/rmats_per_site_psi.csv \
      --splice_table_file output/GRCh38_v29_splice_table.txt \
      --output_file output/rmats/rmats_GRCh38_v29_splice_table.txt
  ```
    - `--rmats_file`: Path to the rMATS per site PSI file.
    - `--splice_table_file`: Path to the splice table file.
    - `--output_file`: Path to the output file.
- Output:
    - A splice table with PSI values and positions from rMATS added to the table.


### SpliSER

`run_spliser.py` 

- Description: Runs SpliSER using the BAM files provided from the bam files TSV and the GTF file. Spliser has three steps
process, combine, and output functions. This script runs all three steps, the combine step is parallelized by chromosome but 
is still slow for large datasets. Regtools is used to extract splice junctions from BAM files. SpliSER is then used to quantify SSE (splice-site strength estimate)
for each splice site of the extracted junctions.
- Run:
```
python src/run_spliser.py --bam_file_list test_data/bam_files.tsv --output_dir output --gtf_file GRCh38/gencode.v29.primary_assembly.annotation.gtf.gz --n_threads 20
```
    - `--bam_file_list`: TSV file containing BAM file paths.
    - `--output_dir`: Directory where output files will be saved.
    - `--gtf_file`: Path to the GTF file.
    - `--n_threads`: Number of threads to use.
    - `--skip_processing` : Skip the processing step and use the processed files in the output directory.
- Output:
  - SpliSER output files are saved in `output_dir/spliser` see SpliSER documentation for more information on the output files.
The final output is spliser_outputAll.DiffSpliSER.tsv which contains the SSE values for each splice site in each BAM file.


`spliser_to_splice_table.py`

- Description: For each transcript in the splice table this scripts adds respective SSE values and genomic positions from SpliSER quantification to the last two columns.
- Run:
```python src/spliser_to_splice_table.py \
    --spliser_output output/spliser/spliser_outputAll.DiffSpliSER.tsv \
    --splice_table output/GRCh38_v29_splice_table.txt \
    --output_file output/spliser/spliser_GRCh38_v29_splice_table.txt \
    --n_threads 20
```
    - `--spliser_ouptut`: Path to the SpliSER DiffSpliSER.tsv file.
    - `--splice_table`: Path to the splice table file.
    - `--output_file`: Path to the output file.
    - `--n_threads`: Number of threads to use.
- Output:
  - A splice table with SSE values and positions from SpliSER added to the last two columns of the table.



    

## References
- LeafCutter: https://github.com/davidaknowles/leafcutter
- rMATS: https://github.com/Xinglab/rmats-turbo
- SpliSER: https://github.com/CraigIDent/SpliSER

