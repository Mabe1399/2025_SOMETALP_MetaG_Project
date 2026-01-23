# BINNING SNAKEFILE
# Last Updated: 2026-01-14
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Binning rules usable in many different pipelines
# The workflow start with the assumption that all the data files are in the same directory.

# Steps included:
# 01: Contig indexing
# 02: Backmapping (all vs all)
# 03: Mixte Binning
#   03a: Metabat2
#   03b: Maxbin2
#   03c: Concoct
# 04: Bin Selection (Das Tool)
# 05: Bin QC (CheckM)
# 06: Bin QC filtering


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

#################### 03 MIXTE BINNING #########################


#################### 03A METABAT2 ####################

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
        basename = lambda wc: f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Bins/{wc.contig}_metabat2/{wc.contig}_metabat2.bin"
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


######################### 03B MAXBIN2 ################################

rule maxbin2_coverage_table:
    input:
        bam = f"{SCRATCH_DIR}/Backmapping/{{contig}}/{{reads}}_mapped_to_{{contig}}_contigs.bam"
    output:
        coverage_table = f"{SCRATCH_DIR}/Coverage_table/{{contig}}/{{reads}}_mapped_to_{{contig}}_coverage.txt"
    conda:
        "../envs/bowtie2_mapping.yaml"
    log:
        f"{LOG_DIR}/03_Binning/MaxBin2/Coverage_table/{{reads}}_mapped_to_{{contig}}_coverage.log"
    shell:
        """
        samtools coverage {input.bam} \
        | tail -n +2 \
        | sort -k1 \
        | cut -f1,6 > {output.coverage_table} 2> {log}
        """


rule maxbin2_abund_list:
    input:
        lambda wc: expand(
            f"{SCRATCH_DIR}/Coverage_table/{wc.contig}/{{reads}}_mapped_to_{wc.contig}_coverage.txt",
            reads=CONTIGS
        )
    output:
        abund_list = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/Abundance_list/{{contig}}_abund_list.txt"
    benchmark:
        f"{BENCH_DIR}/03_Binning/MaxBin2/Abundance_list/{{contig}}_abund_list.tsv"
    log:
        f"{LOG_DIR}/03_Binning/MaxBin2/Abundance_list/{{contig}}_abund_list.log"
    run:
        with open(output.abund_list, 'w') as f:
            for fp in input:
                f.write('%s\n' % fp)


rule MaxBin2:
    input:
        contigs = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        ),
        abund_list = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/Abundance_list/{{contig}}_abund_list.txt"
    output:
        bins = directory(f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/Bins/{{contig}}_maxbin2")
    params:
        basename = lambda wc: f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/Bins/{wc.contig}_maxbin2/{wc.contig}_maxbin2.bin",
        min_contig_length = config["MaxBin2"]["min_contig_length"]
    threads: config["MaxBin2"]["threads"]
    conda:
        "../envs/Binning.yaml"
    benchmark:
        f"{BENCH_DIR}/03_Binning/MaxBin2/Binning/{{contig}}_maxbin2.log"
    log:
        f"{LOG_DIR}/03_Binning/MaxBin2/Binning/{{contig}}_maxbin2.log"
    resources:
        mem_mb = config["MaxBin2"]["memory_mb"],
        runtime = config["MaxBin2"]["runtime_min"],
        cpus_per_task = config["MaxBin2"]["threads"]
    shell:
        """
        mkdir -p {output.bins}

        run_MaxBin.pl -thread {threads} \
        -min_contig_length {params.min_contig_length} \
        -contig {input.contigs} \
        -abund_list {input.abund_list} \
        -out {params.basename} 2> {log} 1>&2

        # Normalize MaxBin2 bin names
        for f in {output.bins}/*.fasta; do
            bn=$(basename "$f" .fasta)
            id=$(echo "$bn" | awk -F. '{{printf "%d", $NF}}')
            mv "$f" "{output.bins}/{{wc.contig}}_maxbin2.bin.${{id}}.fa"
        done

        """


######################### 03C CONCOCT ################################

rule cut_up_fasta:
    input:
        contigs = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        )
    output:
        bed = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/contigs_10K/{{contig}}_contigs.bed",
        contigs_10K = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/contigs_10K/{{contig}}_contigs.fa"
    conda:
        "../envs/Binning.yaml"
    params:
        chunk_size = config["CONCOCT"]["chunk_size"],
        overlap_size = config["CONCOCT"]["overlap_size"]
    benchmark:
        f"{BENCH_DIR}/03_Binning/Concoct/cut_up_fasta/{{contig}}_cut_up.tsv"
    log:
        f"{LOG_DIR}/03_Binning/Concoct/cut_up_fasta/{{contig}}_cut_up.log"
    shell:
        """
        cut_up_fasta.py {input.contigs} \
        -c {params.chunk_size} \
        -o {params.overlap_size} \
        --merge_last \
        -b {output.bed} > {output.contigs_10K} 2> {log}
        """


rule Concoct_coverage_table:
    input:
        bed = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/contigs_10K/{{contig}}_contigs.bed",
        bam = lambda wc: expand(
            f"{SCRATCH_DIR}/Backmapping/{wc.contig}/{{reads}}_mapped_to_{wc.contig}_contigs.bam",
            reads=CONTIGS
        ),
        bai = lambda wc: expand(
            f"{SCRATCH_DIR}/Backmapping/{wc.contig}/{{reads}}_mapped_to_{wc.contig}_contigs.bam.bai",
            reads=CONTIGS)
    output:
        coverage_table = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/coverage_tables/{{contig}}_coverage_table.txt"
    conda:
        "../envs/Binning.yaml"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Concoct/coverage_table/{{contig}}_coverage.tsv"
    log:
        f"{LOG_DIR}/03_Binning/Concoct/coverage_table/{{contig}}_coverage.log"
    resources:
        runtime = config["CONCOCT"]["coverage_runtime"]
    shell:
        """
        concoct_coverage_table.py {input.bed} \
        {input.bam} > {output.coverage_table} 2> {log}
        """


rule Concoct:
    input:
        contig_10K = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/contigs_10K/{{contig}}_contigs.fa",
        coverage_table = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/coverage_tables/{{contig}}_coverage_table.txt"
    output:
        clustering = (
            f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/clustering/{{contig}}_concoct/"
            f"{{contig}}_bins_clustering_gt{config["CONCOCT"]["min_contig_length"]}.csv")
    params:
        bins = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/clustering/{{contig}}_concoct/{{contig}}_bins",
        min_contig_length = config["CONCOCT"]["min_contig_length"]
    conda:
        "../envs/Binning.yaml"
    threads: config["CONCOCT"]["threads"]
    log:
        f"{LOG_DIR}/03_Binning/Concoct/clustering/{{contig}}_clustering.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Concoct/clustering/{{contig}}_clustering.tsv"
    resources:
        mem_mb = config["CONCOCT"]["memory_mb"],
        cpus_per_task = config["CONCOCT"]["threads"],
        runtime = config["CONCOCT"]["runtime_min"]
    shell:
        """
        concoct --threads {threads} -l {params.min_contig_length} \
        --composition_file {input.contig_10K} \
        --coverage_file {input.coverage_table} \
        -b {params.bins}
        2> {log} 1>&2
        """


rule merge_cutup_clustering:
    input:
        bins = (
            f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/clustering/{{contig}}_concoct/"
            f"{{contig}}_bins_clustering_gt{config["CONCOCT"]["min_contig_length"]}.csv")
    output:
        merged = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/merge_cutup_clustering/{{contig}}_clustering_merged.csv"
    conda:
        "../envs/Binning.yaml"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Concoct/merge_cutup_clustering/{{contig}}_clustering_merged.tsv"
    log:
        f"{LOG_DIR}/03_Binning/Concoct/merge_cutup_clustering/{{contig}}_clustering_merged.log"
    shell:
        """
        merge_cutup_clustering.py {input.bins} > {output.merged} 2> {log}
        """


rule extract_fasta_bins:
    input:
        contigs = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        ),
        clustering_merged = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/merge_cutup_clustering/{{contig}}_clustering_merged.csv"
    output:
        fasta_bins = directory(f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/Bins/{{contig}}_concoct")
    conda:
        "../envs/Binning.yaml"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Concoct/extract_fasta_bins/{{contig}}_extract_fasta.tsv"
    log:
        f"{LOG_DIR}/03_Binning/Concoct/extract_fasta_bins/{{contig}}_extract_fasta.log"
    shell:
        """
        mkdir -p {output.fasta_bins}
        extract_fasta_bins.py \
        {input.contigs} \
        {input.clustering_merged} \
        --output_path {output.fasta_bins} \
        2> {log}
        """

rule normalize_concoct_bins:
    input:
        bins = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/Bins/{{contig}}_concoct"
    output:
        normalized = directory(
            f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/Bins/{{contig}}_concoct_normalized"
        )
    localrule: True
    conda:
        "../envs/Binning.yaml"
    shell:
        """
        mkdir -p {output.normalized}
        i=0
        for f in {input.bins}/*.fa; do
            mv "$f" "{output.normalized}/{wildcards.contig}_concoct.bin.${{i}}.fa"
            i=$((i+1))
        done
        """

######################### 04 BIN SELECTION ###########################

rule metabat2_Scaffolds2Bin:
    input:
        bins = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Bins/{{contig}}_metabat2/"
    output:
        scaffolds2bin = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv"
    conda:
        "../envs/Binning.yaml"
    benchmark:
        f"{BENCH_DIR}/03_Binning/MetaBat2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv"
    log:
        f"{LOG_DIR}/03_Binning/MetaBat2/scaffolds2bin/{{contig}}_Scaffolds2Bin.log"
    shell:
        """
            Fasta_to_Contig2Bin.sh \
            -i {input.bins} \
            -e fa \
            | awk 'BEGIN{{OFS="\t"}} {{print $1,$NF}}' > {output.scaffolds2bin}
        """


rule maxbin2_Scaffolds2Bin:
    input:
        bins = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/Bins/{{contig}}_maxbin2/"
    output:
        scaffolds2bin = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv"
    conda:
        "../envs/Binning.yaml"
    benchmark:
        f"{BENCH_DIR}/03_Binning/MaxBin2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv"
    log:
        f"{LOG_DIR}/03_Binning/MaxBin2/scaffolds2bin/{{contig}}_Scaffolds2Bin.log"
    shell:
        """
        Fasta_to_Contig2Bin.sh \
        -i {input.bins} \
        -e fasta > {output.scaffolds2bin}
        """


rule concoct_Scaffolds2Bin:
    input:
        bins = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/Bins/{{contig}}_concoct_normalized/"
    output:
        scaffolds2bin = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv"
    conda:
        "../envs/Binning.yaml"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Concoct/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv"
    log:
        f"{LOG_DIR}/03_Binning/Concoct/scaffolds2bin/{{contig}}_Scaffolds2Bin.log"
    shell:
        """
        Fasta_to_Contig2Bin.sh \
        -i {input.bins} \
        -e fa > {output.scaffolds2bin}
        """


rule DAS_Tool:
    input:
        metabat2 = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv",
        maxbin2 = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv",
        concoct = f"{OUTPUT_DIR}/01_Analysis/03_Binning/Concoct/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv",
        contigs = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        )
    output:
        out = f"{OUTPUT_DIR}/01_Analysis/03_Binning/DAS_Tool/{{contig}}_DASTool/{{contig}}_DASTool_summary.tsv",
        dir = directory(f"{OUTPUT_DIR}/01_Analysis/03_Binning/DAS_Tool/{{contig}}_DASTool/{{contig}}_DASTool_bins")
    params:
        basename = lambda wc: f"{OUTPUT_DIR}/01_Analysis/03_Binning/DAS_Tool/{wc.contig}_DASTool/{wc.contig}",
        search_engine = config["DASTool"]["search_engine"]
    conda:
        "../envs/Binning.yaml"
    threads: config["DASTool"]["threads"]
    benchmark:
        f"{BENCH_DIR}/03_Binning/Das_Tool/{{contig}}_binning.tsv"
    log:
        f"{LOG_DIR}/03_Binning/Das_Tool/{{contig}}_binning.log"
    resources:
        mem_mb = config["DASTool"]["memory_mb"],
        cpus_per_task = config["DASTool"]["threads"],
        runtime = config["DASTool"]["runtime_min"]
    shell:
        """
        DAS_Tool \
        --bins {input.metabat2},{input.maxbin2},{input.concoct} \
        --contigs {input.contigs} \
        --outputbasename {params.basename} \
        --labels metabat2,maxbin2,concoct \
        --write_bins \
        --write_bin_evals \
        --threads {threads} \
        --search_engine {params.search_engine}

        # Rename all bins to consistent pattern DAS_Tool as inconsistent bin naming
        i=1
        for f in {output.dir}/*.fa; do
            mv -n "$f" "{output.dir}/{wildcards.contig}_DASTool_bin_${{i}}.fa"
            i=$((i+1))
        done
        """


######################### 05 BIN QC ##################################

rule CheckM_QC:
    input:
        summary = f"{OUTPUT_DIR}/01_Analysis/03_Binning/DAS_Tool/{{contig}}_DASTool/{{contig}}_DASTool_summary.tsv",
        dir = f"{OUTPUT_DIR}/01_Analysis/03_Binning/DAS_Tool/{{contig}}_DASTool/{{contig}}_DASTool_bins"
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

######################### 06 BIN QC FILTERING ##################################

rule Bin_filtering_QC_summary:
    input:
        stats = expand(f"{OUTPUT_DIR}/01_Analysis/03_Binning/CheckM_QC/{{contig}}_checkm_QC/qa_summary.tsv", contig=CONTIGS),
        bins = expand(f"{OUTPUT_DIR}/01_Analysis/03_Binning/DAS_Tool/{{contig}}_DASTool/{{contig}}_DASTool_bins", contig=CONTIGS)
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
        -e {params.contamination} >> {log}
        """
