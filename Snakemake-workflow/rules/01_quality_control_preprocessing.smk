# QUALITY CONTROL AND PREPROCESSING SNAKEFILE
# Last Updated: 2025-11-18
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
# 05: Before Taxonomic profiling (Kraken2)
# 06: Host filtering (BWA)
# 07: Preprocessing Summary

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
    threads: config["Concatenate"]["threads"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/lane_concatenation/{{sample}}_concat.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/lane_concatenation/{{sample}}_concat.tsv"
    resources:
        runtime = config["Concatenate"]["runtime_min"]
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
    threads: config["fastqc"]["threads"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Before_QC/{{sample}}_QC.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Before_QC/{{sample}}_QC.tsv"
    resources:
        mem_mb = config["fastqc"]["memory_mb"],
        runtime = config["fastqc"]["runtime_min"],
        cpus_per_task = config["fastqc"]["threads"]
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
    threads: config["multiqc"]["threads"]
    log: 
        f"{LOG_DIR}/01_QC_Preprocessing/Before_QC/Multi_QC.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Before_QC/Multi_QC.tsv"
    resources:
        mem_mb = config["multiqc"]["memory_mb"],
        runtime = config["multiqc"]["runtime_min"],
        cpus_per_task = config["multiqc"]["threads"]
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
    threads: config["trimmomatic"]["threads"]
    params:
        nextera = config["trimmomatic"]["trimmomatic_adapters"],
        trailing= config["trimmomatic"]["trailing"],
        leading= config["trimmomatic"]["leading"],
        min_length= config["trimmomatic"]["minlen"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Trimming/{{sample}}_trimmomatic.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Trimming/{{sample}}_trimmomatic.tsv"
    conda:
        "../envs/Trimmomatic.yaml"
    resources:
        mem_mb = config["trimmomatic"]["memory_mb"],
        runtime = config["trimmomatic"]["runtime_min"],
        cpus_per_task = config["trimmomatic"]["threads"]
    shell:
        "trimmomatic PE -phred33 -summary {log} -threads {threads} {input.R1} {input.R2} \
        {output.R1_paired} {output.R1_unpaired} {output.R2_paired} {output.R2_unpaired} \
        ILLUMINACLIP:{params.nextera}:2:30:10 \
        LEADING:{params.leading} TRAILING:{params.trailing} MINLEN:{params.min_length}"

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
    threads: config["fastqc"]["threads"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/After_QC/{{sample}}_QC.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/After_QC/{{sample}}_QC.tsv"
    resources:
        mem_mb = config["fastqc"]["memory_mb"],
        runtime = config["fastqc"]["runtime_min"],
        cpus_per_task = config["fastqc"]["threads"]
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
    threads: config["multiqc"]["threads"]
    log: 
        f"{LOG_DIR}/01_QC_Preprocessing/After_QC/Multi_QC.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/After_QC/Multi_QC.tsv"
    resources:
        mem_mb = config["multiqc"]["memory_mb"],
        runtime = config["multiqc"]["runtime_min"],
        cpus_per_task = config["multiqc"]["threads"]
    shell:
        "mkdir -p {OUTPUT_DIR}/02_Results/01_QC_Preprocessing/QC_After ;"
        "(multiqc -o {output.outdir} {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After) 2> {log}"

###################### 05 BEFORE TAX PROFILING ######################################

rule build_Kraken2_db:
    output:
        directory("/work/FAC/FBM/DBC/fmazel/gut_evol_stg/shared/Database/kraken2_db")
    params:
        host_list=lambda wildcards: " ".join(config["kraken2_db"]["refseq"])
    threads: config["kraken2_db"]["threads"]
    conda:
        "../envs/Kraken2.yaml"
    resources:
        mem_mb = config["kraken2_db"]["memory_mb"],
        runtime = config["kraken2_db"]["runtime_min"],
        disk_mb = config["kraken2_db"]["disk_mb"],
        cpus_per_task = config["kraken2_db"]["threads"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Kraken2/kraken2_db.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Kraken2/kraken2_db.tsv"
    shell:
        """
        kraken2-build --download-taxonomy --db {output}
        kraken2-build --download-library bacteria --db {output}
        kraken2-build --download-library human --db {output}
        for f in {params.host_list}; do
            kraken2-build --add-to-library $f --db {output} 
        done
        kraken2-build --build --db {output}
        """


rule run_kraken2:
    input:
        R1=f"{SCRATCH_DIR}/Trimmed/{{sample}}_R1_paired.fastq.gz",
        R2=f"{SCRATCH_DIR}/Trimmed/{{sample}}_R2_paired.fastq.gz",
        db="/work/FAC/FBM/DBC/fmazel/gut_evol_stg/shared/Database/kraken2_db"
    output:
        tab=temp(f"{SCRATCH_DIR}/Kraken2/{{sample}}_kraken2.out"),
        rep=temp(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/Kraken2/{{sample}}_kraken2.report")
    conda:
        "../envs/Kraken2.yaml"
    threads: config["kraken2_run"]["threads"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Kraken2/{{sample}}_kraken2_run.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Kraken2/{{sample}}_kraken2_run.tsv"
    resources:
        mem_mb = config["kraken2_run"]["memory_mb"],
        runtime = config["kraken2_run"]["runtime_min"],
        cpus_per_task = config["kraken2_run"]["threads"]
    shell:
        "(kraken2 --use-names --threads {threads} \
        --db {input.db} \
        --fastq-input --report {output.rep} --gzip-compressed \
        --paired {input.R1} {input.R2} \
        > {output.tab}) 2> {log}"


rule parse_kraken2_report:
    input:
        expand(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/Kraken2/{{sample}}_kraken2.report", sample=SAMPLES)
    output:
        f"{OUTPUT_DIR}/02_Results/01_QC_Preprocessing/Kraken2/Combined_kraken2_report.tsv"
    localrule: True
    threads: 1
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Kraken2/Combined_kraken2_parsing.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Kraken2/Combined_kraken2_parsing.log"
    shell:
        "python scripts/Combine_kraken_reports.py -o {output} {input}"


########################## 06 HOST FILTERING ########################################

rule host_filtering:
    input: 
        R1=f"{SCRATCH_DIR}/Trimmed/{{sample}}_R1_paired.fastq.gz",
        R2=f"{SCRATCH_DIR}/Trimmed/{{sample}}_R2_paired.fastq.gz"
    output:
        unmapped_R1=f"{OUTPUT_DIR}/00_Data/02_Clean_Data/{{sample}}_R1_HF.fastq.gz",
        unmapped_R2=f"{OUTPUT_DIR}/00_Data/02_Clean_Data/{{sample}}_R2_HF.fastq.gz",
        refstats=f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/HF/{{sample}}_refstats.out"
    conda:
        "../envs/bwa_mapping.yaml"
    threads: config["host_filtering"]["threads"]
    params:
        refs_list=lambda wildcardse: ",".join(config["host_filtering"]["Host_refs"]),
        xmx= config["host_filtering"]["java_mem_u"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/HF/{{sample}}_host_filtering.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/HF/{{sample}}_host_filtering.log"
    resources:
        mem_mb = config["host_filtering"]["memory_mb"],
        runtime = config["host_filtering"]["runtime_min"],
        cpus_per_task = config["host_filtering"]["threads"]
    shell:
        "bbsplit.ch in1={input.R1} in2={input.R2} \
        ref={params.refs_list} \
        basename= $TMPDIR/{wildcards.sample}_HF_discarded_%.sam \
        refstats={output.refstats} rebuild=t \
        outu1={output.unmapped_R1} outu2={output.unmapped_R2} nzo=f -Xmx{params.xmx} threads={threads}"

rule parse_HF_refstats:
    input:
        expand(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/HF/{{sample}}_refstats.out", sample=SAMPLES)
    output:
        f"{OUTPUT_DIR}/02_Results/01_QC_Preprocessing/HF/HF_refstats.tsv"
    localrule: True
    threads: 1
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/HF/HF_refstats_parse.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/HF/HF_refstats_parse.log"       
    run:
        import pandas as pd
        import os

        dfs = []
        for f in input:
            # Extract sample name
            sample = os.path.basename(f).replace("_refstats.out", "")
            # load file in panda
            df = pd.read_csv(f, sep="\t", comment="#")
            # Insert a new sample column
            df.insert(0, "sample", sample)
            # Append to preexisitng files
            dfs.append(df)
        
        # Small reformating + save as tsv file
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_df.to_csv(output[0], sep="\t", index=False)

############################ 07 PREPROCESSING SUMMARY ##############################


