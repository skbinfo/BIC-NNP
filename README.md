# BIC-NNP
Genomic Analysis Pipeline
This script automates the process of downloading reference genomes, fetching assembly data, running quality assessment (QUAST), and performing genome annotation (Prokka) for a specified organism. It organizes the output for easy access and ensures reproducibility in genomic data processing.
Overview
The pipeline:

Loads configuration from a YAML file (config.yaml).
Downloads reference genomes and GFF files from provided URLs.
Searches for genome assemblies using the NCBI Entrez API.
Downloads assemblies with ncbi-genome-download.
Runs QUAST for quality assessment and Prokka for annotation.
Organizes results into designated directories.

Requirements
Python Libraries

subprocess: Runs shell commands.
json: Parses configuration files.
os: Interacts with the operating system.
Bio.Entrez: Fetches genomic data from NCBI.

External Tools



Tool
Version
Purpose



Biopython
1.83
Interacts with NCBI Entrez database.


ncbi-genome-download
0.3.3
Downloads genome assemblies from NCBI.


QUAST
5.2.0
Evaluates genome assembly quality.


Prokka
1.14.6
Annotates prokaryotic genomes.


Maker
2.31.9
Annotates eukaryotic genomes.


wget
-
Downloads files from the internet.


gunzip
-
Decompresses .gz files.


Installation
# Biopython
pip3 install biopython

# ncbi-genome-download
pip3 install ncbi-genome-download

# QUAST
pip3 install quast

# Prokka (use a separate Conda environment if conflicts arise)
conda install -c bioconda prokka

# Maker
conda install bioconda::maker


Note: If Prokka conflicts with other tools, install it in a separate Conda environment and update the run_prokka_with_env function with the correct environment name.

Running the Pipeline
python3 annotation_script.py

Input File

Configuration File: config.yaml

Configuration File Example
# Query for Mycobacterium tuberculosis
organism_query: "Mycobacterium[Organism]"
email: "your@email.com"
organism: "Mycobacterium tuberculosis"
location: "India"
reference_genome_url: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz"
gff_url: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz"
group: "bacteria"
organism_type: "PK" # PK for prokaryotes, EK for eukaryotes
output_folder: "Mycobacterium"

Commands Used in the Script

wget
wget -P {output_dir1} {url}


Purpose: Downloads reference genomes or other files to output_dir1.


os.makedirs
os.makedirs(output_dir1, exist_ok=True)


Purpose: Creates directories like genomes, assembly, or quast_output.


ncbi-genome-download
ncbi-genome-download -l all -F fasta -o {output_dir2} bacteria --assembly-accessions {gcf_accession}


Purpose: Downloads genome sequences in FASTA format from NCBI.


find and mv
find {output_dir2}/*/* -type f -name '*.fna.gz' -exec mv {} {output_dir2} \;


Purpose: Moves .fna.gz files to the main output_dir2 directory.


gunzip
gunzip {output_dir2}/*.gz {output_dir1}/*.gz && find {output_dir2} -type d -empty -delete


Purpose: Decompresses .gz files and removes empty directories.


quast.py
quast.py {assembly_path}/{accession}*.fna -r {reference_genome} -g {gff_genome} -o quast_output/{accession}


Purpose: Evaluates genome assembly quality and saves results to quast_output/{accession}.


prokka
prokka assembly/{asm_id}.fna -outdir prokka_output/{asm_id} --prefix {asm_id}


Purpose: Annotates genomes and saves results to prokka_output/{asm_id}.


os.rename
os.rename(old_path, new_path)


Purpose: Renames files in output_dir2 for consistency.



Output

Genomes: Stored in output_dir1 (e.g., genomes).
Assemblies: Stored in output_dir2 (e.g., assembly).
QUAST Results: Stored in output_dir3 (e.g., quast_output).
Prokka Annotations: Stored in prokka_output/{asm_id}.

Notes

Ensure all required tools are installed and accessible in your environment.
The configuration file (config.yaml) must be correctly formatted and present in the working directory.
For large datasets, monitor disk space and network connectivity.

