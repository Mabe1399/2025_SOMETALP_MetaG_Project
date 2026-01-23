# TAXONOMIC PROFILING SNAKEFILE
# Last Updated: 2026-01-20
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Taxonomic profiling rules usable in many different pipelines
# The workflow start with some assumptions: 
# - cleaned fasta files are present in the directory stored in the DATA_DIR Variable. (Check SNAKEFILE)
# - A Kraken2 database is available and its path is specified in the config file.
# - Sample names are defined in the config file under SAMPLES variable.

# Steps included:
# 01: Taxonomic profiling (Kraken2)


##################################################
# RULES
##################################################

################ 01 TAXONOMIC PROFILING ##########

rule kraken2_Taxonomic_Profiling:
    input:
        R1=f"{DATA_DIR}/02_Clean_Data/{{sample}}_R1_cleaned.fastq.gz",
        R2=f"{DATA_DIR}/02_Clean_Data/{{sample}}_R2_cleaned.fastq.gz",
        db= config["kraken2_db"]["DB_dir"]
    output:
        tab=temp(f"{SCRATCH_DIR}/Kraken2/{{sample}}_kraken2.out"),
        rep=temp(f"{OUTPUT_DIR}/01_Analysis/02b_Taxonomic_Profiling/Kraken2/{{sample}}_kraken2.report")
    conda:
        "../envs/Kraken2.yaml"
    threads: config["kraken2_run"]["threads"]
    log:
        f"{LOG_DIR}/02b_Taxonomic_Profiling/Kraken2/{{sample}}_kraken2_run.log"
    benchmark:
        f"{BENCH_DIR}/02b_Taxonomic_Profiling/Kraken2/{{sample}}_kraken2_run.tsv"
    resources:
        mem_mb = config["kraken2_run"]["memory_mb"],
        runtime = config["kraken2_run"]["runtime_min"],
        cpus_per_task = config["kraken2_run"]["threads"]
    shell:
        "(kraken2 --use-names --threads {threads} \
        --db {input.db} \
        --report {output.rep} --gzip-compressed \
        --report-zero-counts \
        --paired {input.R1} {input.R2} \
        --output {output.tab}) 2> {log}"


rule parse_kraken2_taxonomic_report:
    input:
        expand(f"{OUTPUT_DIR}/01_Analysis/02b_Taxonomic_Profiling/Kraken2/{{sample}}_kraken2.report", sample=SAMPLES)
    output:
        f"{OUTPUT_DIR}/02_Results/02b_Taxonomic_Profiling/Kraken2/Combined_kraken2_tax_profiling_report.tsv"
    localrule: True
    conda:
        "../envs/Parsing.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/02b_Taxonomic_Profiling/Kraken2/Combined_kraken2_parsing.log"
    benchmark:
        f"{BENCH_DIR}/02b_Taxonomic_Profiling/Kraken2/Combined_kraken2_parsing.tsv"
    shell:
        "(python3 scripts/Combine_kraken_reports.py -o {output} -i {input}) 2> {log}"
