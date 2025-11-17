# QUALITY CONTROL AND PREPROCESSING SNAKEFILE
# Last Updated: 2025-11-17
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# QC and Preprocessing rules usable in many different pipelines
# The workflow start with the assumption that all the data files are in the same directory 
# and the sequencing occured in many lanes.

# Steps included:
# 01: Concatenation (When multiple lane sequencing)
# 02: Before QC (FastQC + MultiQC)
# 03: Trimming (Trimmomatic)
# 04: After QC (FastQC + MultiQC)
# 05: Before Taxonomic profiling (Kraken 2)
# 06: Host filtering (BWA)
# 07: After Taxonomic profiling (Kraken 2)

##################################################
# RULES
##################################################

############ 01 Concat Lanes #####################

rule concat_lanes:
    input:
        R1s=lambda wildcards: sorted(glob.glob(f"{DATA_DIR}/{wildcards.sample}_L*_R1.fastq.gz")),
        R2s=lambda wildcards: sorted(glob.glob(f"{DATA_DIR}/{wildcards.sample}_L*_R2.fastq.gz")),
    output:
        R1=f"{SCRATCH_DIR}/Concatenated/{{sample}}_R1_concat.fastq.gz",
        R2=f"{SCRATCH_DIR}/Concatenated/{{sample}}_R2_concat.fastq.gz",
    threads: 2
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/lane_concatenation/{{sample}}_concat.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/lane_concatenation/{{sample}}_concat.tsv"
    resources:
        runtime = 5
    shell:
        "cat {input.R1s} > {output.R1}; cat {input.R2s} > {output.R2} 2> {log}"

################# 02 BEFORE QC #######################

rule before_fastqc:
    input:
        lambda wildcards: [
            f"{SCRATCH_DIR}/Concatenated/{wildcards.sample}_R1_concat.fastq.gz",
            f"{SCRATCH_DIR}/Concatenated/{wildcards.sample}_R2_concat.fastq.gz"
        ],
    output:
        html1=f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before/{{sample}}_R1_concat_fastqc.html",
        html2=f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before/{{sample}}_R2_concat_fastqc.html",
        zip1=f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before/{{sample}}_R1_concat_fastqc.zip",
        zip2=f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before/{{sample}}_R2_concat_fastqc.zip"
    conda:
        "../envs/QC_env.yaml"
    threads: 2
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Before_QC/{{sample}}_QC.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Before_QC/{{sample}}_QC.tsv"
    resources:
        mem_mb = 1000,
        runtime = 35,
        cpus_per_task = 2
    shell:
        "mkdir -p {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before ;"
        "fastqc -t {threads} -o {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before {input};"

rule before_multiqc:
    input:
        expand(
            f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before/{{sample}}_R{{read}}_concat_fastqc.html", 
            sample=SAMPLES, 
            read=[1,2]
            )
    output:
        outdir = directory(f"{OUTPUT_DIR}/02_Results/01_QC_Preprocessing/QC_Before")
    conda:
        "../envs/QC_env.yaml"
    threads: 1
    log: 
        f"{LOG_DIR}/01_QC_Preprocessing/Before_QC/Multi_QC.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Before_QC/Multi_QC.tsv"
    resources:
        mem_mb = 600,
        runtime = 10,
        cpus_per_task = 1
    shell:
        "mkdir -p {OUTPUT_DIR}/02_Results/01_QC_Preprocessing/QC_Before ;"
        "multiqc -o {output.outdir} {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before"

############### 03 TRIMMING #############################

rule Raw_Trimming:
    input:
        R1=f"{SCRATCH_DIR}/Concatenated/{{sample}}_R1_concat.fastq.gz",
        R2=f"{SCRATCH_DIR}/Concatenated/{{sample}}_R2_concat.fastq.gz",
    output:
        R1_paired=f"{SCRATCH_DIR}/Trimmed/{{sample}}_R1_paired.fastq.gz",
        R2_paired=f"{SCRATCH_DIR}/Trimmed/{{sample}}_R2_paired.fastq.gz",
        R1_unpaired=f"{SCRATCH_DIR}/Trimmed/{{sample}}_R1_unpaired.fastq.gz",
        R2_unpaired=f"{SCRATCH_DIR}/Trimmed/{{sample}}_R2_unpaired.fastq.gz"
    threads: 8
    params:
        nextera = config["trimmomatic_adapters"],
        q= 28,
        min_length= 60
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Trimming/{{sample}}_trimmomatic.log"
    conda:
        "../envs/Trimmomatic.yaml"
    resources:
        mem_mb = 6000,
        runtime = 180,
        cpus_per_task = 8
    shell:
        "trimmomatic PE -phred33 -trimlog {log} -threads {threads} {input.R1} {input.R2} \
        {output.R1_paired} {output.R1_unpaired} {output.R2_paired} {output.R2_unpaired} \
        ILLUMINACLIP:{params.nextera}:2:30:10 \
        LEADING:{params.q} TRAILING:{params.q} MINLEN:{params.min_length}"

############### 04 AFTER QC #####################################

rule After_fastqc:
    input:
        lambda wildcards: [
            f"{SCRATCH_DIR}/Trimmed/{wildcards.sample}_R1_paired.fastq.gz",
            f"{SCRATCH_DIR}/Trimmed/{wildcards.sample}_R2_paired.fastq.gz"
        ],
    output:
        html1=f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After/{{sample}}_R1_paired_fastqc.html",
        html2=f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After/{{sample}}_R2_paired_fastqc.html",
        zip1=f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After/{{sample}}_R1_paired_fastqc.zip",
        zip2=f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After/{{sample}}_R2_paired_fastqc.zip"
    conda:
        "../envs/QC_env.yaml"
    threads: 2
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/After_QC/{{sample}}_QC.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/After_QC/{{sample}}_QC.tsv"
    resources:
        mem_mb = 1000,
        runtime = 35,
        cpus_per_task = 2
    shell:
        "mkdir -p {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After ;"
        "(fastqc -t {threads} -o {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After {input}) 2> {log};"
  

rule After_multiqc:
    input:
        expand(
            f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After/{{sample}}_R{{read}}_paired_fastqc.html", 
            sample=SAMPLES, 
            read=[1,2]
            )
    output:
        outdir = directory(f"{OUTPUT_DIR}/02_Results/01_QC_Preprocessing/QC_After")
    conda:
        "../envs/QC_env.yaml"
    threads: 1
    log: 
        f"{LOG_DIR}/01_QC_Preprocessing/After_QC/Multi_QC.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/After_QC/Multi_QC.tsv"
    resources:
        mem_mb = 600,
        runtime = 10,
        cpus_per_task = 1
    shell:
        "mkdir -p {OUTPUT_DIR}/02_Results/01_QC_Preprocessing/QC_After ;"
        "(multiqc -o {output.outdir} {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After) 2> {log}"

###################### 05 HOST FILTERING ######################################
