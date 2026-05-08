# BACTERIAL PHYLOGENY INFERENCE + CO-DIVERSIFICATION ANALYSIS SNAKEFILE
# Last Updated: 2026-05-06
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Bacterial phylogeny inference + Co-Diversification analysis using the MAGs information

# Steps included:
# 01: Bacterial Markers identification (GTDB-tk identify)
# 02: Bacterial Markers alignement and concatenation (GTDB-tk align)
# 03: MAGs Phylogenetic tree inference (IQTree2)
# 04: Co-Diversification Testing (R script)


##################################################
# RULES
##################################################


################### 01 MAGs BACTERIAL MARKERS IDENTIFICATION ##################

rule Mags_GTDBtk_identify:
    input:
        db= config["GTDB_Tk"]["db"],
        Mag= f"{OUTPUT_DIR}/02_Results/03_Binning/HighQC_Bins/"
    output:
        out= directory(f"{OUTPUT_DIR}/01_Analysis/05b_Co_diversification/gtdb_identify/")
    log: 
        f"{LOG_DIR}/05b_Co_diversification/Mags_GTDBtk_identify/Mags_GTDBtk_identify.log"
    benchmark:
        f"{BENCH_DIR}/05b_Co_diversification/Mags_GTDBtk_identify/Mags_GTDBtk_identify.tsv"
    conda:
        "../envs/Reporting.yaml"
    threads: config["GTDB_Tk"]["threads"]
    resources:
        mem_mb = config["GTDB_Tk"]["memory_mb"],
        runtime = config["GTDB_Tk"]["runtime_min"],
        cpus_per_task = config["GTDB_Tk"]["threads"]
    shell:
        """
        export GTDBTK_DATA_PATH={input.db}
        echo "GTDB_Tk DB is in: $GTDBTK_DATA_PATH" >> {log}

        gtdbtk identify --genome_dir {input.Mag} \
            --out_dir {output.out} \
            --extension fa \
            --prefix HighQC_MAGs \
            --cpus {threads} > {log} 2>&1
        """

################### 02 MAGs BACTERIAL MARKER ALIGNEMENT AND CONCATENATION ##################

rule Mags_GTDBtk_align:
    input:
        db= config["GTDB_Tk"]["db"],
        identify= f"{OUTPUT_DIR}/01_Analysis/05b_Co_diversification/gtdb_identify/"
    output: 
        out= directory(f"{OUTPUT_DIR}/01_Analysis/05b_Co_diversification/gtdb_align/"),
        msa= f"{OUTPUT_DIR}/01_Analysis/05b_Co_diversification/gtdb_align/HighQC_MAGs.bac120.user_msa.fasta.gz"
    log:
        f"{LOG_DIR}/05b_Co_diversification/Mags_GTDBtk_align/Mags_GTDBtk_align.log"
    benchmark:
        f"{BENCH_DIR}/05b_Co_diversification/Mags_GTDBtk_align/Mags_GTDBtk_align.tsv"
    conda:
        "../envs/Reporting.yaml"
    threads: config["GTDB_Tk"]["threads"]
    resources:
        mem_mb = config["GTDB_Tk"]["memory_mb"],
        runtime = config["GTDB_Tk"]["runtime_min"],
        cpus_per_task = config["GTDB_Tk"]["threads"]
    shell:
        """
        export GTDBTK_DATA_PATH={input.db}
        echo "GTDB_Tk DB is in: $GTDBTK_DATA_PATH" >> {log}

        gtdbtk align --identify_dir {input.identify} \
            --out_dir {output.out} \
            --skip_gtdb_refs \
            --prefix HighQC_MAGs \
            --cpus {threads} > {log} 2>&1
        """

####################### 03 MAGs PHYLOGENY INFERENCE ##################

rule iqtree_phylo_infer:
    input:
        msa= f"{OUTPUT_DIR}/01_Analysis/05b_Co_diversification/gtdb_align/HighQC_MAGs.bac120.user_msa.fasta.gz"
    output:
        tree= f"{OUTPUT_DIR}/01_Analysis/05b_Co_diversification/iqtree/HighQC_MAGs.bac120.user_msa.treefile"
    log:
        f"{LOG_DIR}/05b_Co_diversification/MAG_Phylo_inference/iqtree.log"
    benchmark:
        f"{BENCH_DIR}/05b_Co_diversification/MAG_Phylo_inference/iqtree.tsv"
    conda:
        "../envs/Phylogeny.yaml"
    threads: config["Phylo_iqtree"]["threads"]
    resources:
        mem_mb = config["Phylo_iqtree"]["memory_mb"],
        runtime = config["Phylo_iqtree"]["runtime_min"],
        cpus_per_task = config["Phylo_iqtree"]["threads"]
    shell:
        """
        outdir=$(dirname {output.tree})

        iqtree -s {input.msa} \
        --prefix ${{outdir}}/HighQC_MAGs.bac120.user_msa \
        -B 1000 \
        -nt {threads} \
        -T {threads} > {log} 2>&1
        """

