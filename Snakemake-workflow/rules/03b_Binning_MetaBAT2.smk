# BINNING METABAT2 SNAKEFILE
# Last Updated: 2026-01-20
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Metabat2 Binning rules usable in many different pipelines
# The workflow start with the assumption that all the data files are in the same directory.

# Steps included:
# 01: Contig indexing
# 02: Backmapping (all vs all)
# 03: Metabat2
# 04: Bin QC (CheckM)
# 05: Bin QC filtering


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
        bam = temp(f"{SCRATCH_DIR}/Backmapping/{{contig}}/{{reads}}_mapped_to_{{contig}}_contigs.bam"),
        bai = temp(f"{SCRATCH_DIR}/Backmapping/{{contig}}/{{reads}}_mapped_to_{{contig}}_contigs.bam.bai")
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
            2>> {log} \
            | samtools view -@ {threads} -m 1G -bS - 2>> {log} \
            | samtools sort -@ {threads} -m 1G -o {output.bam} 2>> {log}
        
        samtools index {output.bam} 2>> {log}
        """

#################### 03 METABAT2 ####################

rule Depth_calculations:
    input:
       bam = lambda wc: expand(
            f"{SCRATCH_DIR}/Backmapping/{wc.contig}/{{reads}}_mapped_to_{wc.contig}_contigs.bam",
            reads=CONTIGS
        )
    output:
        depth = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Depth_Contig/{{contig}}_depth.txt"
    conda:
        "../envs/Binning.yaml"
    log:
        f"{LOG_DIR}/03_Binning/MetaBat2/Depth/{{contig}}_depth.log"
    shell:
        """
        echo "Using BAMs:" >&2
        printf "  %s\n" {input.bam} >&2

        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth} \
            {input.bam} 2> {log}
        """


rule Metabat2:
    input:
        contig = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        ),
        depth = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Depth_Contig/{{contig}}_depth.txt"
    output:
        Bins = directory(f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Bins/{{contig}}_metabat2")
    params:
        min_contig_length = config["MetaBat2"]["min_contig_length"],
        basename = lambda wc: f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Bins/{wc.contig}_metabat2/{wc.contig}_bin"
    conda:
        "../envs/Binning.yaml"
    threads: config["MetaBat2"]["threads"]
    log:
        f"{LOG_DIR}/03_Binning/MetaBat2/Binning/{{contig}}_metabat2.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/MetaBat2/Binning/{{contig}}_metabat2.tsv"
    resources:
        mem_mb = config["MetaBat2"]["memory_mb"],
        runtime = config["MetaBat2"]["runtime_min"],
        cpus_per_task = config["MetaBat2"]["threads"]
    shell:
        """
        mkdir -p {output.Bins}
        metabat2 --numThreads {threads} \
                 --inFile {input.contig} \
                 --outFile {params.basename} \
                 --abdFile {input.depth} \
                 --minContig {params.min_contig_length} 2> {log} 1>&2
        """

######################### 04 BIN QC ##################################

rule CheckM_QC:
    input:
        dir = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Bins/{{contig}}_metabat2/"
    output:
        file = f"{OUTPUT_DIR}/01_Analysis/03_Binning/CheckM_QC/{{contig}}_checkm_QC/qa_summary.tsv"
    params:
        outdir = lambda wc: f"{OUTPUT_DIR}/01_Analysis/03_Binning/CheckM_QC/{wc.contig}_checkm_QC",
        db = config["CheckM"]["db"]
    conda:
        "../envs/Binning.yaml"
    log:
        f"{LOG_DIR}/03_Binning/CheckM/{{contig}}_checkm_QC.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/CheckM/{{contig}}_checkm_QC.tsv"
    threads: config["CheckM"]["threads"]
    resources:
        mem_mb = config["CheckM"]["memory_mb"],
        runtime = config["CheckM"]["runtime_min"],
        cpus_per_task = config["CheckM"]["threads"]
    shell:
        """
        export CheckM_DATA_PATH={params.db}
        checkm lineage_wf {input.dir} {params.outdir} -x fa -t {threads}
        checkm qa {params.outdir}/lineage.ms {params.outdir} -f {output.file} -t {threads} -o 1 --tab_table
        """

######################### 05 BIN QC FILTERING ##################################

rule Bin_filtering_QC_summary:
    input:
        stats = expand(f"{OUTPUT_DIR}/01_Analysis/03_Binning/CheckM_QC/{{contig}}_checkm_QC/qa_summary.tsv", contig=CONTIGS),
        bins = expand(f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Bins/{{contig}}_metabat2", contig=CONTIGS)
    output:
        full_stats= f"{OUTPUT_DIR}/02_Results/03_Binning/CheckM_QC/all_MAGs_stats.tsv",
        filtered_stats= f"{OUTPUT_DIR}/02_Results/03_Binning/CheckM_QC/filtered_MAGs_stats.tsv",
        filtered_mags= directory(f"{OUTPUT_DIR}/02_Results/03_Binning/HighQC_Bins/")
    conda:
        "../envs/Parsing.yaml"
    log:
        f"{LOG_DIR}/03_Binning/CheckM/aggregate_checkm.log"
    threads: 1
    localrule: True
    params:
        completeness= config["Bin_filtering"]["completeness"],
        contamination= config["Bin_filtering"]["contamination"]
    shell:
        """
        python scripts/aggregate_checkm.py -i {input.stats} \
        -b {input.bins} \
        -o {output.full_stats} \
        -f {output.filtered_stats} \
        -m {output.filtered_mags} \
        -c {params.completeness} \
        -e {params.contamination}
        """
