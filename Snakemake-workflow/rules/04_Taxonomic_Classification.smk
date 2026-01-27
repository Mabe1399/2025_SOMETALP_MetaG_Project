# TAXONOMIC CLASSIFICATION SNAKEFILE
# Last Updated: 2026-01-20
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Taxonomic Classification rules usable in many different pipelines to classify MAGs created

# Steps included:
# 01: Taxonomic Classification (GTDB-Tk)
# 02: iTOL Parsing


##################################################
# RULES
##################################################

########## 01 TAXONOMIC CLASSIFICATION ###########

rule gtdb_classify:
    input:
        db= config["GTDB_Tk"]["db"],
        filtered_mags= f"{OUTPUT_DIR}/02_Results/03_Binning/HighQC_Bins/"
    output:
        classified_out= directory(f"{OUTPUT_DIR}/01_Analysis/04_Taxonomic_Classification/gtdb_classify/")
    log:
        f"{LOG_DIR}/04_Taxonomic_Classification/GTDB_Tk/GTDB_classify.log"
    benchmark:
        f"{BENCH_DIR}/04_Taxonomic_Classification/GTDB_Tk/GTDB_classify.tsv"
    conda:
        "../envs/Taxonomy.yaml"
    threads: config["GTDB_Tk"]["threads"]
    params:
        tmp = f"{SCRATCH_DIR}/gtdbtk_classify/",
        scratch = f"{SCRATCH_DIR}/gtdbtk_classify_test_pplacer"
    resources:
        mem_mb= config["GTDB_Tk"]["memory_mb"],
        cpus_per_task= config["GTDB_Tk"]["threads"],
        runtime= config["GTDB_Tk"]["runtime_min"]
    shell:
        """
        export GTDBTK_DATA_PATH={input.db}
        echo "GTDB_Tk DB is in: $GTDBTK_DATA_PATH" >> {log}

        mkdir -p {params.tmp}
        
        gtdbtk classify_wf --genome_dir {input.filtered_mags} \
        --extension fa \
        --tmpdir {params.tmp} \
        --pplacer_cpus 4 \
        --scratch_dir {params.scratch} \
        --out_dir {output.classified_out} \
        --cpus {threads} 2> {log} 1>&2
        """

################ 02 iTOL PARSING #########################

rule convert_to_iTOL:
    input:
        db = config["GTDB_Tk"]["db"],
        dir = f"{OUTPUT_DIR}/01_Analysis/04_Taxonomic_Classification/gtdb_classify/",
        tree = f"{OUTPUT_DIR}/01_Analysis/04_Taxonomic_Classification/gtdb_classify/classify/gtdbtk.bac120.classify.tree.{{cls}}.tree"
    output:
        iTOL = f"{OUTPUT_DIR}/02_Results/04_Taxonomic_Classification/iTOL_visualisation/gtdbtk.bac120.iTOL.{{cls}}.tree"
    log:
        f"{LOG_DIR}/04_Taxonomic_Classification/iTOL_visualisation/iTOL_tree_{{cls}}_formating.log"
    conda:
        "../envs/Taxonomy.yaml"
    threads: 1
    localrule: True
    shell:
        """
        # Still need the DB
        export GTDBTK_DATA_PATH={input.db}
        echo "GTDB_Tk DB is in: $GTDBTK_DATA_PATH" >> {log}

        # Reformat the tree to be usable in iTOL
        gtdbtk convert_to_itol --input_tree {input.tree} --output_tree {output.iTOL} 2> {log} 1>&2
        
        """

rule parse_gtdbtk_summary_to_iTOL:
    input:
        gtdb_summary = f"{OUTPUT_DIR}/01_Analysis/04_Taxonomic_Classification/gtdb_classify/classify/gtdbtk.bac120.summary.tsv"
    output:
        dir = directory(f"{OUTPUT_DIR}/02_Results/04_Taxonomic_Classification/iTOL_visualisation/Annotation_files")
    params:
        metadata= config["DATA_VARIABLES"]["metadata_path"],
        config= "../config/config.yaml"
    log:
        f"{LOG_DIR}/04_Taxonomic_Classification/iTOL_visualisation/iTOL_Annotation_files.log"
    conda:
        "../envs/Parsing.yaml"
    threads: 1
    localrule: True
    shell:
        """
        python3 scripts/parse_gtdbtk_to_itol.py --metadata {params.metadata} \
        --summary {input.gtdb_summary} \
        --config {params.config} \
        --outdir {output.dir} 2> {log} 1>&2
        """