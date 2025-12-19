# ASSEMBLY SNAKEFILE
# Last Updated: 2025-12-19
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Assembly rules usable in many different pipelines
# The workflow start with the assumption that all the data files are in the same directory.

# Steps included:
# 01: Assembly (Metaspade)
# 02: Contig QC (Metaquast)
# 03: filter out low quality contig
# 04: Assembly Summary

##################################################
# RULES
##################################################

############## 01 ASSEMBLY #######################

rule Metaspade_assembly:
    input:
        R1 = f"{OUTPUT_DIR}/00_Data/02_Clean_Data/{{sample}}_R1_cleaned.fastq.gz",
        R2 = f"{OUTPUT_DIR}/00_Data/02_Clean_Data/{{sample}}_R2_cleaned.fastq.gz"
    output:
        scaffolds = f"{OUTPUT_DIR}/02_Results/02_Assembly/Metaspade/{{sample}}_contigs.fasta"
    conda:
        "../envs/Assembly.yaml"
    params:
        dir= directory("$TMPDIR/{sample}_metaspades/")
    logs:
        f"{LOG_DIR}/02_Assembly/Metaspade/{{sample}}_contigs.log"
    benchmark:
        f"{BENCH_DIR}/02_Assembly/Metaspade/{{sample}}_contigs.tsv"
    threads: config["Metaspade_assembly"]["threads"]
    resources:
        mem_mb = config["Metaspade_assembly"]["memory_mb"],
        runtime = config["Metaspade_assembly"]["runtime_min"],
        cpus_per_task = config["Metaspade_assembly"]["threads"]
    shell:
        """
        spades.py --meta \
            --pe-1 {input.R1} --pe1-2 {input.R2} \
            -o {params.dir} \
            -m $(({resources.mem_mb}/1024)) \
            -t {threads} 2> {log}
        
        # Move the file of interest to the right directory
        mv {params.dir}/contig.fasta {output.scaffolds}
        """


################ 02 CONTIG QC #############################

rule Metaquast_QC:
    input:
        f"{OUTPUT_DIR}/02_Results/02_Assembly/Metaspade/{{sample}}_contigs.fasta"
    output:
        directory(f"{OUTPUT_DIR}/01_Analysis/02_Assembly/Metaquast_QC/{{sample}}")
    conda:
        "../envs/Assembly.yaml"
    logs:
        f"{LOG_DIR}/02_Assembly/MetaQuast/{{sample}}_QC.log"
    benchmark:
        f"{BENCH_DIR}/02_Assembly/MetaQuast/{{sample}}_QC.tsv"
    threads: config["MetaQuast"]["threads"]
    resources:
        mem_mb = config["MetaQuast"]["memory_mb"],
        runtime = config["MetaQuast"]["runtime_min"],
        cpus_per_task = config["MetaQuast"]["threads"]
    shell:
        """
        metaquast.py \
            {input} \
            -t {threads} \
            -o {output}
        """

rule Quast_multiqc:
    input:
        expand(f"{OUTPUT_DIR}/01_Analysis/02_Assembly/Metaquast_QC/{{sample}}", sample=SAMPLES)
    output:
        f"{OUTPUT_DIR}/02_Results/02_Assembly/Metaquast_QC/multiqc_report.html"
    conda:
        "../envs/QC_env.yaml"
    log: 
        f"{LOG_DIR}/02_Assembly/MetaQuast/Multi_QC.log"
    benchmark:
        f"{BENCH_DIR}/02_Assembly/MetaQuast/Multi_QC.tsv"
    threads: config["multiqc"]["threads"]
    resources:
        mem_mb = config["multiqc"]["memory_mb"],
        runtime = config["multiqc"]["runtime_min"],
        cpus_per_task = config["multiqc"]["threads"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/02_Results/02_Assembly/Metaquast_QC
        multiqc {OUTPUT_DIR}/01_Analysis/02_Assembly/Metaquast_QC \
            --outdir {OUTPUT_DIR}/02_Results/02_Assembly/Metaquast_QC \
            --force
        """
    
################ 03 FILTER CONTIG ############################