# BINNING SNAKEFILE
# Last Updated: 2026-01-05
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Binning rules usable in many different pipelines
# The workflow start with the assumption that all the data files are in the same directory.

# Steps included:
# 01: Contig indexing
# 02: Backmapping (all vs all)
# 03: Depth calculations (Metabat2)
# 04: Binning (Metabat2)
# 05: Bin QC (CheckM)


##################################################
# RULES
##################################################

#################### 01 CONTIG INDEXING ####################

rule Contig_indexing:
    input:
        contig = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        )
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
        prefix = lambda wc: f"{SCRATCH_DIR}/Contig_Index/{wc.contig}/{wc.contig}_index",
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
        bowtie2-build --threads {threads} {input.contig} {params.prefix} 2> {log}
        """

#################### 02 BACKMAPPING ####################

rule Contig_Backmapping:
    input:
        R1 = f"{OUTPUT_DIR}/00_Data/02_Clean_Data/{{reads}}_R1_cleaned.fastq.gz",
        R2 = f"{OUTPUT_DIR}/00_Data/02_Clean_Data/{{reads}}_R2_cleaned.fastq.gz",
        index = [
            f"{SCRATCH_DIR}/Contig_Index/{{contig}}/{{contig}}_index.1.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{contig}}/{{contig}}_index.2.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{contig}}/{{contig}}_index.3.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{contig}}/{{contig}}_index.4.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{contig}}/{{contig}}_index.rev.1.bt2",
            f"{SCRATCH_DIR}/Contig_Index/{{contig}}/{{contig}}_index.rev.2.bt2",
        ]
    output:
        bam = temp(f"{SCRATCH_DIR}/Backmapping/{{contig}}/{{reads}}_mapped_to_{{contig}}_contigs.bam")
    params:
        index = lambda wc: f"{SCRATCH_DIR}/Contig_Index/{wc.contig}/{wc.contig}_index"
    conda:
        "../envs/bowtie2_mapping.yaml"
    threads: config["Contig_Backmapping"]["threads"]
    log:
        f"{LOG_DIR}/03_Binning/Mapping/{{reads}}_to_{{contig}}_Contig_backmapping.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Mapping/{{reads}}_to_{{contig}}_Contig_backmapping.tsv"
    resources:
        mem_mb = config["Contig_Backmapping"]["memory_mb"],
        runtime = config["Contig_Backmapping"]["runtime_min"],
        cpus_per_task = config["Contig_Backmapping"]["threads"]
    shell:
        """
        mkdir -p $(dirname {output.bam})
        bowtie2 --threads {threads} \
            --no-unal \
            -x {params.index} \
            -1 {input.R1} \
            -2 {input.R2} \
            | samtools view -@ {threads} -m 1G -bS - \
            | samtools sort -@ {threads} -m 1G -o {output.bam} 2> {log}
        
        samtools index {output.bam}
        """


#################### 03 DEPTH CALCULATIONS ####################

rule Depth_calculations:
    input:
       bam = lambda wc: expand(
            f"{SCRATCH_DIR}/Backmapping/{wc.contig}/{{reads}}_mapped_to_{wc.contig}_contigs.bam",
            reads=SAMPLES
        )
    output:
        depth = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Depth_Contig/{{contig}}_depth.txt"
    conda:
        "../envs/Binning.yaml"
    threads: config["Depth_calc"]["threads"]
    log:
        f"{LOG_DIR}/03_Binning/Depth/{{contig}}_depth.log"
    shell:
        """
        echo "Using BAMs:" >&2
        printf "  %s\n" {input.bam} >&2

        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth} \
            {input.bam} 2> {log}
        """


#################### 04 BINNING ####################

rule Metabat2:
    input:
        contig = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        ),
        depth = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Depth_Contig/{{contig}}_depth.txt"
    output:
        dir = directory(f"{OUTPUT_DIR}/02_Results/03_Binning/Bins/{{contig}}_metabat2/")
    params:
        min_contig_length = config["Metabat2"]["min_contig_length"],
        basename = lambda wc: f"{OUTPUT_DIR}/02_Results/03_Binning/Bins/{wc.contig}_metabat2/{wc.contig}_bin"
    conda:
        "../envs/Binning.yaml"
    threads: config["Metabat2"]["threads"]
    log:
        f"{LOG_DIR}/03_Binning/Binning/{{contig}}_binning.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Binning/{{contig}}_binning.tsv"
    resources:
        mem_mb = config["Metabat2"]["memory_mb"],
        runtime = config["Metabat2"]["runtime_min"],
        cpus_per_task = config["Metabat2"]["threads"]
    shell:
        """
        mkdir -p {output.dir}
        metabat2 --numThreads {threads} \
                 --inFile {input.contig} \
                 --outFile {params.basename} \
                 --abdFile {input.depth} \
                 --minContig {params.min_contig_length} 2> {log} 1>&2
        """

######################### 05 BIN QC ##################################

