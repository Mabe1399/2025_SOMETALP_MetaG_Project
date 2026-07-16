#!/usr/bin/env python3
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)
# 
# Usage: python3 parse_MAG_representative_info.py --CheckM_summary CheckM_summary_file.tsv \
# --GTDB_Tk_summary GTDB_Tk_summary_file.tsv \
# --drep_dir drep_directory \
# --output output_file.tsv

# Output file: MAG_representative_info_parsed.tsv
# This script parses the CheckM summary, GTDB-Tk summary and dRep directory to extract information about the representative MAGs and others, including their genome ID, completeness, contamination, taxonomy and dRep cluster.
# The output file will contain the following columns: MAG ID, sample, completeness, contamination, Genome size (bp), N50, # markers, dRep_cluster, taxonomy. 

# IMPORTS
import pandas as pd
import argparse
import os

# PARSE ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument("--CheckM_summary", required=True)
parser.add_argument("--GTDB_Tk_summary", required=True)
parser.add_argument("--drep_dir", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

# LOAD CheckM SUMMARY
checkm_summary = pd.read_csv(args.CheckM_summary, sep="\t")

# LOAD GTDB-Tk SUMMARY
gtdbtk_summary = pd.read_csv(args.GTDB_Tk_summary, sep="\t", na_values="N/A")

# EXTRACT GENOME INFO FROM CheckM SUMMARY
## Extract genome ID, completeness and contamination
genome_info = checkm_summary[["Bin Id", "Completeness", "Contamination", "# markers" ]]
genome_info.columns = ["MAG ID", "completeness", "contamination", "# markers"]
# Add .fasta extension to genome ID
genome_info.loc[:, "MAG ID"] = genome_info["MAG ID"].apply(lambda x: x + ".fa")

# EXTRACT TAXONOMY INFO FROM GTDB-Tk SUMMARY
## Extract genome ID and taxonomy
taxonomy_info = gtdbtk_summary[["user_genome", "classification"]]
taxonomy_info.columns = ["MAG ID", "taxonomy"]
# Add .fasta extension to genome ID
taxonomy_info.loc[:, "MAG ID"] = taxonomy_info["MAG ID"].apply(lambda x: x + ".fa")
# widen taxonomy column to separate different taxonomic ranks
taxonomy_info[["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]] = taxonomy_info["taxonomy"].str.split(";", expand=True)
# drop original taxonomy column
taxonomy_info = taxonomy_info.drop(columns=["taxonomy"])

# EXTRACT dRep CLUSTER INFO FROM dRep DIRECTORY
## Extract info from data_tables directory
## Extract Score from Sdb.csv file
drep_data_tables_dir = os.path.join(args.drep_dir, "data_tables")
sdb_file = os.path.join(drep_data_tables_dir, "Sdb.csv")
sdb_df = pd.read_csv(sdb_file)
drep_cluster_info = sdb_df[["genome", "score"]]
drep_cluster_info.columns = ["MAG ID", "dRep_Score"]

## Extract cluster info from Cdb.csv file
cdb_file = os.path.join(drep_data_tables_dir, "Cdb.csv")
cdb_df = pd.read_csv(cdb_file)
drep_cluster_info_cdb = cdb_df[["genome", "secondary_cluster"]]
drep_cluster_info_cdb.columns = ["MAG ID", "dRep_cluster"]

## Extract genome size, N50, closest cluster member, furthest cluster member from Widb.csv file
widb_file = os.path.join(drep_data_tables_dir, "Widb.csv")
widb_df = pd.read_csv(widb_file)
drep_cluster_info_widb = widb_df[["genome","closest_cluster_member", "furthest_cluster_member"]]
drep_cluster_info_widb.columns = ["MAG ID","closest_cluster_member", "furthest_cluster_member"]

# MERGE ALL DATAFRAMES if info NoT available for a MAG, fill with NA
## Merge genome_info with taxonomy_info
merged_df = pd.merge(genome_info, taxonomy_info, on="MAG ID", how="left")
## Merge with drep_cluster_info
merged_df = pd.merge(merged_df, drep_cluster_info, on="MAG ID", how="left")
## Merge with drep_cluster_info_cdb
merged_df = pd.merge(merged_df, drep_cluster_info_cdb, on="MAG ID", how="left")
## Merge with drep_cluster_info_widb
merged_df = pd.merge(merged_df, drep_cluster_info_widb, on="MAG ID", how="left")

# SAVE OUTPUT
merged_df.to_csv(args.output, sep="\t", index=False)