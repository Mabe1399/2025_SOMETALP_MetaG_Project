# QUALITY CONTROL AND PREPROCESSING SNAKEFILE
# Last Updated: 2025-11-10
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
    resources:
        runtime_s = 600
    shell:
        "(cat {input.R1s} > {output.R1}; cat {input.R2s} > {output.R2}) 2> {log}"

################# BEFORE QC #######################

rule before_fastqc:
    input:
        expand(f"{SCRATCH_DIR}/Concatenated/{{sample}}_R{{read}}_concat.fastq.gz", sample=SAMPLES, read=["1", "2"])
    output:
        expand(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before/{{sample}}_R{{read}}_concat_fastqc.html", sample=SAMPLES, read=["1","2"]),
        expand(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before/{{sample}}_R{{read}}_concat_fastqc.zip", sample=SAMPLES, read=["1","2"])
    threads: 4
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Before_QC.log"
    resources:
        mem_mb = 6000,
        runtime= 3600
    shell:
        "mkdir -p {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before; "
        "fastqc -t {threads} -o {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before {input};"