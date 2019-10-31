# Isoseq_processing
Snakemake for aligning Isoseq reads and processing transcripts using Cupcake TOFU and SQANTI

For use in linux

Developed with snakemake (v3.13.3)

# Create the conda environment 
conda env create -f environment_sqanti.yaml

# Install auxillary programs and create directory structure
Install Cupcake/ToFU (https://github.com/Magdoll/cDNA_Cupcake/wiki)
Note: This requires a python3 environment

Install SQANTI (https://bitbucket.org/ConesaLab/sqanti/src/master/), not SQANTI2
(https://github.com/Magdoll/SQANTI2)
Note: Both SQANTI and SQANTI2 require a python2 environment, but have different dependencies.

Data files (e.g. polished isoseq transcripts and cluster reports) can be 
located within a different directory from the directory used to hold results.
This functionality is intended to prevent unnecessary file duplications.

Currently, the Snakemake file will combine the input fasta files and process the
single file. Further modifications are underway to allow individual fasta files
to be processed independently.

Running the Snakemake process on one or more Isoseq fasta files produces the results files from 
ToFU and creates a shell script for running SQANTI. Further modification of the Snakefile is needed
to allow SQANTI to successfully be called from the Snakefile.

```
.
├── config.yaml
├── cluster.json
├── run_snakemake_cluster.sh
└── Snakefile
```

To run the snakefile:

conda activate py3 ##activate python3 environment containing snakemake

snakemake -nr --snakefile Snakefile ##shows what commands will be run

nohup sh run_snakemake_cluster.sh & ## to run commands via bsub cluster system
