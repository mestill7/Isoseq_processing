configfile: "config.yaml" ##Add config["merge_samples"] == "TRUE"

SAMPLES = open(config["sample_list"]) ##
SAMPLES = SAMPLES.read()
SAMPLES = SAMPLES.split() ## a list

import os

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

createFolder(expand("{outdir}/align",outdir=config["out_dir"])[0])
createFolder(expand("{outdir}/temp",outdir=config["out_dir"])[0])

rule all:
    input:
        expand("{outdir}/align/Merged_samples.sorted.sam",outdir=config["out_dir"])[0],
        expand("{outdir}/align/Merged_samples_gmap.collapsed.gff",outdir=config["out_dir"])[0],
        expand("{outdir}/align/Merged_samples_gmap.collapsed.abundance.txt",outdir=config["out_dir"])[0],
        expand("{outdir}/temp/generate_sqanti_qc.sh",outdir=config["out_dir"])[0]

rule merge_samples:
    input:
        expand("{fqdir}/{sample}.fasta",sample=SAMPLES,fqdir=config["fastq_dir"])
    output:
        expand("{outdir}/Merged_samples.fasta",outdir=config["out_dir"])[0]
    run:
        shell("cat {input} > {output}")

rule gmap_align:
    input:
        fa = expand("{outdir}/Merged_samples.fasta",outdir=config["out_dir"])[0],
        gindex_location = expand("{ginloc}",ginloc=config["genome_index_loc"])[0]
    output:
        temp(expand("{outdir}/align/Merged_samples.sam",outdir=config["out_dir"])[0])
    params:
        config["genome_index_name"]
    log:
        "logs/Merged_samples.alignment.log"
    shell:
        "module load gmap/2019-02-15\n"
        "gmap -D {input.gindex_location} -d {params} -f samse -n 0 -t 15 --max-intronlength-ends 200000 -z sense_force {input.fa} > {output}"

rule sort_gmap:
    input:
        expand("{outdir}/align/Merged_samples.sam",outdir=config["out_dir"])[0]
    output:
        expand("{outdir}/align/Merged_samples.sorted.sam",outdir=config["out_dir"])[0]
    params:
        gname = config["genome_index_name"],
        out_dir = config["out_dir"]
    run:
        shell("samtools sort -T {params.out_dir}/temp/{params.gname} -o {output} {input}")

rule collapse_sam:
    input:
        fasta = expand("{outdir}/Merged_samples.fasta",outdir=config["out_dir"])[0],
        sam = expand("{outdir}/align/Merged_samples.sorted.sam",outdir=config["out_dir"])[0]
    output:
        expand("{outdir}/align/Merged_samples_gmap.collapsed.gff",outdir=config["out_dir"])[0]
    params:
        config["out_dir"]
    run:
        shell("module load gmap/2019-02-15\n"),
        shell("module load R\n"),
        shell("PATH=$PATH:/sc/hydra/work/estilm01/Tools/cDNA_Cupcake/sequence/\n"),
        shell("which collapse_isoforms_by_sam.py\n"),
        shell("collapse_isoforms_by_sam.py --input {input.fasta} -s {input.sam} --dun-merge-5-shorter -o {params}/align/Merged_samples_gmap")

rule get_abundance:
    input:
        cluster = expand("{clusterrep}",clusterrep=config["cluster_report"])[0],
        collapsed = expand("{outdir}/align/Merged_samples_gmap.collapsed.gff",outdir=config["out_dir"])[0]
    output:
        expand("{outdir}/align/Merged_samples_gmap.collapsed.abundance.txt",outdir=config["out_dir"])[0]
    params:
        config["out_dir"]
    run:
        shell("module load gmap/2019-02-15\n"),
        shell("module load R\n"),
        shell("PATH=$PATH:/sc/hydra/work/estilm01/Tools/cDNA_Cupcake/sequence/\n"),
        shell("get_abundance_post_collapse.py {params}/align/Merged_samples_gmap.collapsed {input.cluster}")

rule generate_squanti_script:
    input:
        collapsed_gff = expand("{outdir}/align/Merged_samples_gmap.collapsed.gff",outdir=config["out_dir"])[0],
        fl = expand("{outdir}/align/Merged_samples_gmap.collapsed.abundance.txt",outdir=config["out_dir"])[0]
    params:
        index = config["genome_index_loc"],
        index_name = config["genome_index_name"],
        ref_gtf = config["gtf"],
        genome_seq = config["genome_sequence"],
        out_dir = config["out_dir"]
    output:
        expand("{outdir}/temp/generate_sqanti_qc.sh",outdir=config["out_dir"])[0]
    run:
        shell("echo \"conda activate sqanti\nmodule load gmap/2019-02-15\nmodule load R\nexport PYTHONPATH=\$PYTHONPATH:/path/to/Tools/cDNA_Cupcake/sequence/\" > {output}"),
        shell("echo \"cd {params.out_dir}/align\" >> {output}"),
        shell("echo \"python /path/to/Tools/sqanti/sqanti_qc.py -g -x {params.index}/{params.index_name} -t 20 -o Merged_samples_sqanti --fl_count {input.fl} {input.collapsed_gff} {params.ref_gtf} {params.genome_seq}\" >> {output}"),
        shell("echo \"cd {params.out_dir}/align\" >> {output}"),
        shell("echo \"python /path/to/Tools/sqanti/sqanti_filter.py {params.out_dir}/align/Merged_samples_sqanti_classification.txt\" >> {output}"),
        shell("echo \"conda deactivate\" >> {output}")
