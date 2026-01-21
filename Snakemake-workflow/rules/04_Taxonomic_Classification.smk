# TAXONOMIC CLASSIFICATION SNAKEFILE
# Last Updated: 2026-01-20
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Taxonomic Classification rules usable in many different pipelines to classify MAGs created

# Steps included:
# 01: Taxonomic Classification (GTDB-Tk)
# 02: 


##################################################
# RULES
##################################################

########## 01 TAXONOMIC CLASSIFICATION ###########

rule gtdb_db_deploy:
    output:
        dir= directory(config["GTDB_Tk"]["db"])
    log:
        f"{LOG_DIR}/04_Taxonomic_Classification/GTDB_Tk/database_deploy.log"
    localrule: True
    shell:
        """
        echo "GTDB database directory creation" >> {log}
        mkdir -p {output.dir}
        echo "Get Archived Database" >> {log}
        wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
        echo "Unarchive the DB" >> {log}
        tar xvzf gtdbtk_data.tar.gz
        echo "Database deployed and ready to use" >> {log}
        """


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
    ressources:
        mem_mb= config["GTDB_Tk"]["memory_mb"],
        cpus_per_task= config["GTDB_Tk"]["threads"],
        runtime= config["GTDB_Tk"]["runtime_min"]
    shell:
        """
        export GTDBTK_DATA_PATH={input.db}
        echo "GTDB_Tk DB is in: $GTDBTK_DATA_PATH" >> {log}

        gtdbtk classify_wf --genome_dir {input.filtered_mags} \
        --skip_ani_screen \
        --extension fa \
        --out_dir {output.classified_out} \
        --cpus {threads} 2> {log} 1>&2
        """