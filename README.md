# Isoseq_processing
Snakemake for aligning Isoseq reads and processing transcripts using Cupcake TOFU and SQANTI

For use in linux

Developed with snakemake (v3.13.3)

# Create the conda environment 
conda env create -f environment_sqanti.yaml

Data files (e.g. polished isoseq transcripts and cluster reports) can be 
located within a different directory from the directory used to hold results.
This functionality is intended to prevent unnecessary file duplications.

Currently, the Snakemake file will combine the input fasta files and process the
single file. Further modifications are underway to allow individual fasta files
to be processed independently.


# Install auxillary programs and create directory structure

Install ngsutilsj
(https://github.com/compgen-io/ngsutilsj)

Note: rule sort_combined_lengths will require a correct path for ngsutilsj

conda activate py3 ##activate python3 environment containing snakemake

snakemake -nr --snakefile Snakefile ##shows what commands will be run

nohup sh run_snakemake_cluster.sh & ## to run commands via bsub cluster system
