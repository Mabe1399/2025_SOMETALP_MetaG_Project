# INSTRAIN MICROBIAL DIVERSITY ANALYSIS SNAKEFILE
# Last Updated: 2026-02-25
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Instrain microbial diversity analysis rules usable for different pipelines.

# Steps included:
# 01: SGBs Catalogue Annotation
# 02: SGBs Catalogue Mapping (Bowtie2)
# 03: Sample Profiling (Instrain)
# 04: Sample Comparison (Instrain)


##################################################
# RULES
##################################################

############## 01 SGBs CATALOGUE ANNOTATION ###############

rule SGB_catalogue_file:
    input:
        dRep= f"{OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/reference_SGBs/"
    output:
        rep= temp(f"{OUTPUT_DIR}/01_Analysis/05a_Instrain_Analysis/SGBs_representatives_genomes.fasta")
    log:
        f"{LOG_DIR}/05a_Instrain_Analysis/SGBs_representatives_genomes.log"
    conda:
        "../envs/Parsing.yaml"
    threads: 1
    localrule: True
    shell:
        """
        cat {input.dRep}/dereplicated_genomes/* > {output.rep}
        """

rule Scaffold_to_bin_info:
    input:
        dRep= f"{OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/reference_SGBs/"
    output:
        rep= temp(f"{OUTPUT_DIR}/01_Analysis/05a_Instrain_Analysis/SGBs_representatives_genomes.stb")
    log:
        f"{LOG_DIR}/05a_Instrain_Analysis/SGBs_Scaffolds_to_bin.log"
    conda:
        "../envs/Reporting.yaml"
    threads: 1
    localrule: True
    shell:
        """
        parse_stb.py --reverse -f {input.dRep}/dereplicated_genomes/* -o {output.rep} 2> {log} 1>&2
        """

rule prodigal_SGBs_prediction:
    input:
        rep= f"{OUTPUT_DIR}/01_Analysis/05a_Instrain_Analysis/SGBs_representatives_genomes.fasta"
    output:
        faa= temp(f"{OUTPUT_DIR}/01_Analysis/05a_Instrain_Analysis/SGBs_representatives_genomes.genes.faa"),
        fna= temp(f"{OUTPUT_DIR}/01_Analysis/05a_Instrain_Analysis/SGBs_representatives_genomes.genes.fna")
    log:
        f"{LOG_DIR}/05a_Instrain_Analysis/SGBs_gene_prediction_prodigal.log"
    conda:
        "../envs/Annotation.yaml"
    threads: 1
    localrule: True
    shell:
        """
        prodigal -i {input.rep} \
        -d {output.fna} \
        -a {output.faa} -p meta 2> {log} 1>&2
        """


#################### 02 SGBs CATALOGUE MAPPING ##########################

rule SGBs_representatives_indexing:
    input:
        SGBs = f"{OUTPUT_DIR}/01_Analysis/05a_Instrain_Analysis/SGBs_representatives_genomes.fasta"
    output:
        temp([
            f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index.1.bt2",
            f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index.2.bt2",
            f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index.3.bt2",
            f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index.4.bt2",
            f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index.rev.1.bt2",
            f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index.rev.2.bt2",
        ])
    params:
        prefix = f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index",
        outdir = f"{SCRATCH_DIR}/SGBs_Catalogue_index"
    conda:
        "../envs/bowtie2_mapping.yaml"
    threads: config["SGBs_indexing"]["threads"]
    log:
        f"{LOG_DIR}/05a_Instrain_Analysis/SGBs_representatives_genomes_indexing.log"
    benchmark:
        f"{BENCH_DIR}/05a_Instrain_Analysis/SGBs_representatives_genomes_indexing.tsv"
    resources:
        mem_mb = config["SGBs_indexing"]["memory_mb"],
        runtime = config["SGBs_indexing"]["runtime_min"],
        cpus_per_task = config["SGBs_indexing"]["threads"]
    shell:
        """
        mkdir -p {params.outdir}
        bowtie2-build --threads {threads} {input.SGBs} {params.prefix} 2> {log} 1>&2
        """

rule SGBs_Sample_mapping:
    input:
        R1 = f"{OUTPUT_DIR}/00_Data/02_Clean_Data/{{reads}}_R1_cleaned.fastq.gz",
        R2 = f"{OUTPUT_DIR}/00_Data/02_Clean_Data/{{reads}}_R2_cleaned.fastq.gz",
        index = [
            f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index.1.bt2",
            f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index.2.bt2",
            f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index.3.bt2",
            f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index.4.bt2",
            f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index.rev.1.bt2",
            f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index.rev.2.bt2",
        ]
    output:
        bam = temp(f"{SCRATCH_DIR}/SGBs_mapping/{{reads}}/{{reads}}_mapped_to_SGBs_representatives.bam"),
        bai = temp(f"{SCRATCH_DIR}/SGBs_mapping/{{reads}}/{{reads}}_mapped_to_SGBs_representatives.bam.bai")
    params:
        index = f"{SCRATCH_DIR}/SGBs_Catalogue_index/SGBs_representatives_genomes_index"
    conda:
        "../envs/bowtie2_mapping.yaml"
    threads: config["SGBs_mapping"]["threads"]
    log:
        f"{LOG_DIR}/05a_Instrain_Analysis/SGBs_mapping/{{reads}}_to_SGBs_mapping.log"
    benchmark:
        f"{BENCH_DIR}/05a_Instrain_Analysis/SGBs_mapping/{{reads}}_to_SGBs_mapping.tsv"
    resources:
        mem_mb = config["SGBs_mapping"]["memory_mb"],
        runtime = config["SGBs_mapping"]["runtime_min"],
        cpus_per_task = config["SGBs_mapping"]["threads"]
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

################# 03 SAMPLE PROFILING #####################

rule Sample_Profiling_instrain:
    input:
        bam= f"{SCRATCH_DIR}/SGBs_mapping/{{reads}}/{{reads}}_mapped_to_SGBs_representatives.bam",
        rep= f"{OUTPUT_DIR}/01_Analysis/05a_Instrain_Analysis/SGBs_representatives_genomes.fasta",
        stb= f"{OUTPUT_DIR}/01_Analysis/05a_Instrain_Analysis/SGBs_representatives_genomes.stb",
        fna= f"{OUTPUT_DIR}/01_Analysis/05a_Instrain_Analysis/SGBs_representatives_genomes.genes.fna"
    output:
        profile= directory(f"{OUTPUT_DIR}/02_Results/05a_Instrain_Analysis/Instrain_profile/{{reads}}_mapped_to_SGBs_representatives.IS")
    conda:
        "../envs/Instrain.yaml"
    threads: config["InStrain"]["threads"]
    log:
        f"{LOG_DIR}/05a_Instrain_Analysis/Instrain_profile/{{reads}}_mapped_to_SGBs_representatives.log"
    benchmark:
        f"{BENCH_DIR}/05a_Instrain_Analysis/Instrain_profile/{{reads}}_mapped_to_SGBs_representatives.tsv"
    resources:
        mem_mb = config["InStrain"]["memory_mb"],
        runtime = config["InStrain"]["runtime_min"],
        cpus_per_task = config["InStrain"]["threads"]
    shell:
        """
        inStrain profile {input.bam} {input.rep} \
        -p {threads}  \
        -g {input.fna} \
        -s {input.stb} \
        -o {output.profile} \
        --database_mode 2> {log} 1>&2
        """

rule Sample_Compare_instrain:
    input:
        profile= expand(f"{OUTPUT_DIR}/02_Results/05a_Instrain_Analysis/Instrain_profile/{{reads}}_mapped_to_SGBs_representatives.IS/", reads = CONTIGS),
        stb= f"{OUTPUT_DIR}/01_Analysis/05a_Instrain_Analysis/SGBs_representatives_genomes.stb"
    output:
        compare= directory(f"{OUTPUT_DIR}/02_Results/05a_Instrain_Analysis/Instrain_compare/SGBs_representatives_genomes.COMPARE")
    conda:
        "../envs/Instrain.yaml"
    threads: config["InStrain"]["threads"]
    log:
        f"{LOG_DIR}/05a_Instrain_Analysis/Instrain_compare/SGBs_representatives_genomes.COMPARE.log"
    benchmark:
        f"{BENCH_DIR}/05a_Instrain_Analysis/Instrain_compare/SGBs_representatives_genomes.COMPARE.tsv"
    resources:
        mem_mb = config["InStrain"]["memory_mb"],
        runtime = config["InStrain"]["runtime_min"],
        cpus_per_task = config["InStrain"]["threads"]
    shell:
        """
        inStrain compare \
        -i {input.profile} \
        -p {threads} \
        -s {input.stb} \
        -o {output.compare} \
        --database_mode 2> {log} 1>&2
        """

rule parse_Instrain_output:
    input:
        Instrain_dir= expand(f"{OUTPUT_DIR}/02_Results/05a_Instrain_Analysis/Instrain_profile/{{reads}}_mapped_to_SGBs_representatives.IS", reads= CONTIGS),
        stb = f"{OUTPUT_DIR}/01_Analysis/05a_Instrain_Analysis/SGBs_representatives_genomes.stb",
        dRep = f"{OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/MAG_Representative_info_parsed.tsv"
    output:
        out = directory(f"{OUTPUT_DIR}/02_Results/05a_Instrain_Analysis/Instrain_profile/Parsed")
    log:
        f"{LOG_DIR}/05a_Instrain_Analysis/Instrain_profile/parse_Instrain_profile.log"
    conda:
        "../envs/Parsing.yaml"
    threads: 1
    localrule: True
    shell:
        """
        python3 scripts/parse_Instrain_profile.py --instrain_output {input.Instrain_dir} \
        --MAG_representative_info {input.dRep} \
        --stb {input.stb} \
        --output_dir {output.out} 2> {log} 1>&2
        """

