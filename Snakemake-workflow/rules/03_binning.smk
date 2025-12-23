# BINNING SNAKEFILE
# Last Updated: 2025-12-23
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Binning rules usable in many different pipelines
# The workflow start with the assumption that all the data files are in the same directory.

# Steps included:
# 01: Backmapping (all vs all)
# 02: Depth calculations (Metabat2)
# 03: Binning (Metabat2)
# 04: Bin QC (CheckM)


##################################################
# RULES
##################################################

################# 01 BACKMAPPING #################

rule Contig_indexing:
    input:
        f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/{{contig}}_min{{minlen}}_contigs.fasta"
    output:
        temp([
            f"{SCRATCH_DIR}/Contig_Index/{{contig}}/{{contig}}_index.1.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{contig}}/{{contig}}_index.2.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{contig}}/{{contig}}_index.3.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{contig}}/{{contig}}_index.4.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{contig}}/{{contig}}_index.rev.1.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{contig}}/{{contig}}_index.rev.2.bt2",
        ])
    params:
        prefix = lambda wc: f"{SCRATCH_DIR}/Contig_Index/{wc.contig}/{wc.contig}_index"
        outdir = lambda wc: f"{SCRATCH_DIR}/Contig_Index/{wc.contig}"
    conda:
        "../envs/bowtie2_mapping.yaml"
    threads: config["Contig_indexing"]["threads"]
    log:
        f"{LOG_DIR}/03_Binning/Indexing/{{contig}}_Contig_indexing.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Indexing/{{contig}}_Contig_indexing.tsv"
    resources:
        mem_mb = config["Contig_indexing"]["memory_mb"],
        runtime = config["Contig_indexing"]["runtime_min"],
        cpus_per_task = config["Contig_indexing"]["threads"]
    shell:
        """
        mkdir -p {params.outdir}
        # Build the index
        (bowtie2-build --threads {threads} {input} {params.prefix}) 2> {log}
        """

rule Contig_Backmapping:
    input:
        R1 = lambda wc: f"{OUTPUT_DIR}/00_Data/02_Clean_Data/{{wc.sample}}_R1_cleaned.fastq.gz",
        R2 = lambda wc: f"{OUTPUT_DIR}/00_Data/02_Clean_Data/{{wc.sample}}_R2_cleaned.fastq.gz",
        index = lambda wc: [
            f"{SCRATCH_DIR}/Contig_Index/{{wc.contig}}/{{wc.contig}}_index.1.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{wc.contig}}/{{wc.contig}}_index.2.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{wc.contig}}/{{wc.contig}}_index.3.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{wc.contig}}/{{wc.contig}}_index.4.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{wc.contig}}/{{wc.contig}}_index.rev.1.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{wc.contig}}/{{wc.contig}}_index.rev.2.bt2",
        ]
    output:
        bam=temp(f"{SCRATCH_DIR}/Backmapping/{{wc.contig}}/{{wc.sample}}_mapped_to_{{wc.contig}}_contigs.bam")
    conda:
        "../envs/bowtie2_mapping.yaml"
    threads: config["Contig_Backmapping"]["threads"]
    params:
        index = lambda wc: f"{SCRATCH_DIR}/Contig_Index/{wc.contig}/{wc.contig}_index"
    log:
        f"{LOG_DIR}/03_Binning/Mapping/{{sample}}_to_{{contig}}_Contig_backmapping.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Mapping/{{sample}}_to_{{contig}}_Contig_backmapping.tsv"
    resources:
        mem_mb = config["Contig_Backmapping"]["memory_mb"],
        runtime = config["Contig_Backmapping"]["runtime_min"],
        cpus_per_task = config["Contig_Backmapping"]["threads"]
    shell:
        """
        mkdir -p $(dirname {output.bam})
        (bowtie2 --threads {threads} \
            --no-unal \
            -x {params.index} \
            -1 {input.R1} \
            -2 {input.R2} \
        | samtools view -@ {threads} -b \
        | samtools sort -@ {threads} -o {output.bam}) 2> {log}
        """


################### 02 DEPTH CALCULATIONS #############################

rule Depth_calculations:
    input:
        bam = lambda wc: expand(
            f"{SCRATCH_DIR}/Backmapping/{{wc.contig}}/{{sample}}_mapped_to_{{wc.contig}}_contigs.bam", 
            sample=SAMPLES
        )
    output:
        depth = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Depth_Contig/{{contig}}_depth.txt"
    conda:
        "../envs/Binning.yaml"
    threads: config["Depth_calc"]["threads"]
    log:
        f"{LOG_DIR}/03_Binning/Depth/{{contig}}_depth.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Depth/{{contig}}_depth.tsv"
    resources:
        mem_mb = config["Depth_calc"]["memory_mb"],
        runtime = config["Depth_calc"]["runtime_min"],
        cpus_per_task = config["Depth_calc"]["threads"]
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam} 2> {log}
        """