# Isoseq_processing
Snakemake for aligning Isoseq reads and processing transcripts using Cupcake TOFU and SQANTI

For use in linux

Developed with snakemake (v3.13.3)

# Create the conda environment 
conda env create -f environment_isoseq.yaml

### Test which version of R you are working with!
R (>= 3.4.0)
R packages (for sqanti_qc2.py): ggplot2, scales, reshape, gridExtra, grid, dplyr

# Install auxillary programs and create directory structure
### Install Cupcake/ToFU (https://github.com/Magdoll/cDNA_Cupcake/wiki)
Basic steps:

conda activate isoseq_proc

cd ~/Tools/

git clone https://github.com/Magdoll/cDNA_Cupcake.git

cd ~/Tools/cDNA_Cupcake

python setup.py build

python setup.py install

### Install TAMA

cd ~/Tools/

git clone https://github.com/mestill7/tama.git

### Install SQANTI2
https://github.com/Magdoll/SQANTI2/

Basic steps:

cd \~/Tools

git clone https://github.com/Magdoll/SQANTI2.git


### Install auxillary script (gtfToGenePred)
cd \~/Tools

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred


### Install Matchannot
https://github.com/TomSkelly/MatchAnnot.git

Basic steps:

cd \~/Tools

git clone https://github.com/TomSkelly/MatchAnnot.git


## Once all prerequisites are installed:
### Initialize the environment paths
conda activate isoseq_proc

export PYTHONPATH=$PYTHONPATH:~/Tools/cDNA_Cupcake/sequence/

export PYTHONPATH=$PYTHONPATH:~/Tools/cDNA_Cupcake/

module load R

export PATH=$PATH:\~/Tools


### The following commands should now work (error-free):
python ~/Tools/SQANTI2/sqanti_qc2.py --help

python ~/Tools/SQANTI2/sqanti_filter2.py --help

### Create directory structure for snakemake process
For this example, let's say your directory is "\~/isoseq"

Copy the following files into your "\~/isoseq" directory,

such that your final directory structure is as follows:

```
.
├── config.yaml
├── cluster.json
├── run_snakemake_cluster.sh
├── run_snakemake.sh
└── Snakefile
```
Please modify the config.yaml as needed. Note that: 
Data files (e.g. polished isoseq transcripts and cluster reports) can be 
located within a different directory from the directory used to hold results.
This functionality is intended to prevent unnecessary file duplications.

# Running the Snakefile:
_Make sure to initialize your environment paths as shown previously_
snakemake -nr --snakefile Snakefile ##shows what commands will be run

#### To use your LSF clusters:
```nohup sh run_snakemake_cluster.sh &```

#### To run Snakemake without an LSF cluster:
Type ```sh run_snakemake.sh``` followed by the maximum number of CPU cores to be
used by snakemake. For example, type ```sh run_snakemake.sh 2``` for 2 CPU cores. 
You can also type ```nohup sh run_snakemake.sh 2 &``` to run the pipeline in background.

Currently, the Snakemake file will combine the input fasta files and process the
single file. Future implementations of the pipeline will allow individual fasta 
files to be processed independently.

Running the Snakemake process produces the results files from TAMA and 
runs SQANTI2 and matchannot.  

### Simplifying matchannot results
The matchannot program produces an extremely detailed report. However, users may
only be interested in the associated MatchAnnot score for a given Isoseq isoform.
In that case, the following awk script will extract the isoform name (1st column), 
associated known transcript (2nd, 3rd columns) and MatchAnnot score (4th column):

_Enter directory holding your matchannot results first_

```awk 'BEGIN {OFS="\t"}; {if($1=="result:") print $2,$3,$4,$8}' Merged_samples_matchannot.txt > Merged_samples_matchannot_simplified.txt```
