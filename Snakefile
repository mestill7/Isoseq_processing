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
        expand("{outdir}/align/Merged_samples_sqanti2_completion",outdir=config["out_dir"])[0]

rule merge_samples:
    input:
        expand("{fqdir}/{sample}.fasta",sample=SAMPLES,fqdir=config["fastq_dir"])
    output:
        expand("{outdir}/Merged_samples.fasta",outdir=config["out_dir"])[0]
    run:
        shell("cat {input} > {output}")

rule gmap_build_idx:
    input:
        expand("{ginloc}",ginloc=config["genome_sequence"])[0]
    output:
        expand("{ginloc}/{gname}/{gname}.chromosome",ginloc=config["genome_index_loc"],gname = config["genome_index_name"])[0]
    log:
        "logs/Create_gmap_index.log"
    params:
        gindex_location = config["genome_index_loc"],
        gname = config["genome_index_name"]
    shell:
        "gmap_build -d {params.gname} -D {params.gindex_location} {input}"

rule gmap_align:
    input:
        fa = expand("{outdir}/Merged_samples.fasta",outdir=config["out_dir"])[0],
        gindex_veritas = expand("{ginloc}/{gname}/{gname}.chromosome",ginloc=config["genome_index_loc"],gname = config["genome_index_name"])[0]
    output:
        temp(expand("{outdir}/align/Merged_samples.sam",outdir=config["out_dir"])[0])
    params:
        gloc = config["genome_index_loc"],
        gname = config["genome_index_name"]
    log:
        "logs/Merged_samples.alignment.log"
    shell:
        "gmap -D {params.gloc} -d {params.gname} -f samse -n 0 -t 15 --max-intronlength-ends 200000 -z sense_force {input.fa} > {output}"

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
        bed1 = expand("{outdir}/align/Merged_samples_tama.bed",outdir=config["out_dir"])[0],
        bed2 = expand("{outdir}/align/Merged_samples_tama_trans_read.bed",outdir=config["out_dir"])[0]
    params:
        outdir = config["out_dir"],
        genome_seq = config["genome_sequence"],
        tama_dir = config["tama_dir"]
    run:
        shell("{params.tama_dir}/tama_collapse.py -s {input.sam} -f {params.genome_seq} -p {params.outdir}/align/Merged_samples_tama -x no_cap -i 95 -z 1000")    
        # shell("python {params.tama_dir}/tama_collapse.py -s {input.sam} -f {params.genome_seq} -p {params.outdir}/align/Merged_samples_tama -x no_cap -i 95 -z 1000")
        
# Nancy step 2: change the TAMA collapsed isoform IDs to match SQANTI2’s expected collapsed isoform annotations (PB.XXXX):
rule reformat_bed:
    input:
        bed = expand("{outdir}/align/Merged_samples_tama.bed",outdir=config["out_dir"])[0]
    output:
        expand("{outdir}/align/Merged_samples_tama_corrected.fasta",outdir=config["out_dir"])[0],
        expand("{outdir}/align/Merged_samples_tama_corrected.bed",outdir=config["out_dir"])[0]
    params:
        outdir = config["out_dir"],
        genome_seq = config["genome_sequence"]
    run:
        import pandas as pd
        with open(input.bed, "r") as infile:
            fin_list=list()
            for lines in infile:
                lines = lines.strip().split("\t")
                toadd = [lines[0],lines[1],lines[2],lines[3].split(";")[1].replace('G','PB.')]
                toadd.extend(lines[4:12])
                fin_list.append(toadd)

        df = pd.DataFrame(fin_list)  # df[3].replace('G','PB.',inplace=True)
        df.to_csv("align/Merged_samples_tama_corrected.bed", sep='\t',header=False,index=False)

        shell("bedtools getfasta -name -split -s -fi {params.genome_seq} \
                 -bed align/Merged_samples_tama_corrected.bed -fo align/Merged_samples_tama_corrected.fasta")
        shell("sed -i 's/::.*//g' align/Merged_samples_tama_corrected.fasta") # Format fasta line headers for sqanti
        shell("sed -i 's/^>G/>PB./g' align/Merged_samples_tama_corrected.fasta") # Format fasta line headers for sqanti
        # shell("awk '{gsub(\"G\",\"PB.\",$4)}1' align/Merged_samples_tama_corrected.bed \
        #     > align/Merged_samples_tama_corrected_fixed.bed")
        # shell("mv align/Merged_samples_tama_corrected_fixed.bed align/Merged_samples_tama_corrected.bed")

rule create_group:
    input:
        bed = expand("{outdir}/align/Merged_samples_tama_trans_read.bed",outdir=config["out_dir"])[0]
    output:
        expand("{outdir}/align/Merged_samples.collapsed.group.txt",outdir=config["out_dir"])[0]
    params:
        out_dir = config["out_dir"]
    run:
        import pandas as pd
        with open(input.bed, "r") as infile:
            s_list=list()
            for lines in infile:
                lines = lines.strip().split("\t")
                s_list.append(lines[3].split(";"))
        df = pd.DataFrame(s_list)
        df.columns=['txt','rna']
        df_comp = df.groupby(['txt'], as_index = False).agg({'rna': ','.join}) 
        df_comp['txt'] = [x.replace('G','PB.') for x in df_comp['txt']] 
        df_comp.to_csv(output[0], sep='\t',header=False,index=False)

rule get_abundance:
    input:
        cluster = expand("{clusterrep}",clusterrep=config["cluster_report"])[0],
        collapsed = expand("{outdir}/align/Merged_samples_tama_corrected.fasta",outdir=config["out_dir"])[0],
        group = expand("{outdir}/align/Merged_samples.collapsed.group.txt",outdir=config["out_dir"])[0]
    output:
        expand("{outdir}/align/Merged_samples.collapsed.abundance.txt",outdir=config["out_dir"])[0]
    params:
        config["out_dir"]
    run:
        shell("get_abundance_post_collapse.py {params}/align/Merged_samples.collapsed {input.cluster}")

rule convertbed2gtf:
    input:
        bed = expand("{outdir}/align/Merged_samples_tama_corrected.bed",outdir=config["out_dir"])[0]
    output:
        expand("{outdir}/align/Merged_samples_tama_corrected.gtf",outdir=config["out_dir"])[0]
    params:
        outdir = config["out_dir"],
        genome_seq = config["genome_sequence"]
    run:
        shell("bedToGenePred {input.bed} Merged_samples_tama_genepred")
        shell("genePredToGtf file Merged_samples_tama_genepred {output}")
        shell("rm Merged_samples_tama_genepred")

rule generate_squanti:
    input:
        collapsed_fa = expand("{outdir}/align/Merged_samples_tama_corrected.fasta",outdir=config["out_dir"])[0],
        collapsed_gtf = expand("{outdir}/align/Merged_samples_tama_corrected.gtf",outdir=config["out_dir"])[0],        
        fl = expand("{outdir}/align/Merged_samples.collapsed.abundance.txt",outdir=config["out_dir"])[0]
    params:
        index = config["genome_index_loc"],
        index_name = config["genome_index_name"],
        ref_gtf = config["gtf"],
        genome_seq = config["genome_sequence"],
        out_dir = config["out_dir"],
        sq_dir = config["sqanti_dir"],
        ma_dir = config["matchannot_dir"]
    output:
        expand("{outdir}/align/Merged_samples_sqanti2_completion",outdir=config["out_dir"])[0]
    run:
        shell("python2 {params.ma_dir}/matchAnnot.py --gtf={params.ref_gtf} \
            --format='alt' {params.out_dir}/align/Merged_samples.sorted.sam  > \
            {params.out_dir}/align/Merged_samples_matchannot.txt")
        shell("python {params.sq_dir}/sqanti_qc2.py --gtf -x {params.index}/{params.index_name} -t 30 \
            -d {params.out_dir}/align -o Merged_samples_sqanti --fl_count {input.fl} \
            {input.collapsed_gtf} {params.ref_gtf} {params.genome_seq}")
        shell("echo \"SQANTI2 QC complete\" > {params.out_dir}/align/Merged_samples_sqanti2_completion")
        shell("python {params.sq_dir}/sqanti_filter2.py --skipGTF \
            --faa {params.out_dir}/align/Merged_samples_sqanti_corrected.faa \
            --sam {params.out_dir}/align/Merged_samples.sorted.sam \
            {params.out_dir}/align/Merged_samples_sqanti_classification.txt \
            {input.collapsed_fa} {input.collapsed_gtf}")

