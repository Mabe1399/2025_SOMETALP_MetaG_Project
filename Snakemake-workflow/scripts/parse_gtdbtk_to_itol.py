#!/usr/bin/env python3
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)
# 
# Usage: python parse_gtdbtk_to_itol.py --summary gtdbtk_classify_output --metadata metadata_file --config config.yaml --outdir output_itol_dir

# Output file: leaf_genus_gtdb_to_iTOL_strip_labels.txt
# Output file: leaf_genus_gtdb_to_iTOL_labels.txt
# Output file: leaf_genus_gtdb_to_iTOL_heatmap.txt

# IMPORTS
import pandas as pd
import yaml
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import argparse

# PARSE ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument("--metadata", required=True)
parser.add_argument("--summary", required=True)
parser.add_argument("--config", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

# LOAD METADATA
if args.metadata.endswith(".xlsx"):
    metadata = pd.read_excel(args.metadata)
else:
    metadata = pd.read_csv(args.metadata, sep="\t")

# LOAD GTDB-TK SUMMARY
gtdb_summary = pd.read_csv(args.summary, sep="\t")

# LOAD CONFIG COLORS
with open(args.config) as f:
    cfg = yaml.safe_load(f)

HOST_COLUMN_NAME = cfg["DATA_VARIABLES"]["col_host"]
SPLIT_CRITERIA = cfg["DATA_VARIABLES"]["split_criteria"]
SAMPLE_ID_COL = cfg["DATA_VARIABLES"]["col_sample_id"]

HOST_COLORS = cfg.get("host_colors", None)

#########################################################

# DATA MERGE AND PROCESSING
## Parse user_genome to match metadata
gtdb_summary["Sample_ID"] = gtdb_summary["user_genome"].str.split(SPLIT_CRITERIA).str[0]

## Merge metadata with summary
merged = pd.merge(gtdb_summary, metadata, right_on=SAMPLE_ID_COL, left_on="Sample_ID", how="left")

# LEAF ID EXTRACTION
## Extract Class/ Genus from GTDB-Tk taxonomy
def parse_tax(t):
    tax_levels = t.split(";")
    tax = {}
    for level in tax_levels:
        if "__" in level:
            key, value = level.split("__", 1)
            tax[key] = value
    return tax

## Create a taxonomic table
tax_table = merged["classification"].apply(parse_tax).apply(pd.Series)
## set empty strings to NaN
tax_table = tax_table.replace("", pd.NA)


## extract leaf names from tax_table
merged["leaf_genus"] = "g__" + tax_table["g"].astype(str)
# if genus is NaN, use family
merged["leaf_genus"] = merged["leaf_genus"].fillna("f__" + tax_table["f"].astype(str))

# LEAF LABEL CREATION
## Create leaf ID to leaf label mapping
leaf_id_to_label = dict(zip(merged["leaf_genus"],merged["leaf_genus"]))

## change "f__" labels to "Unclassified_<family_name>"
for leaf_id, label in leaf_id_to_label.items():
    if label.startswith("f__"):
        family_name = label[3:]
        leaf_id_to_label[leaf_id] = f"Unclassified_{family_name}"

## change "g__" labels to just genus name
for leaf_id, label in leaf_id_to_label.items():
    if label.startswith("g__"):
        genus_name = label[3:]
        leaf_id_to_label[leaf_id] = genus_name

## add number of genomes per leaf to label
leaf_counts = merged["leaf_genus"].value_counts()
for leaf_id in leaf_id_to_label.keys():
    count = leaf_counts.get(leaf_id, 0)
    leaf_id_to_label[leaf_id] += f" (n={count})"

# HEATMAP CREATION

## presence/absence of leaf genera across samples in a matrix (0/1  format)
heatmap_data = pd.crosstab(merged[SAMPLE_ID_COL], merged["leaf_genus"])
heatmap_data = (heatmap_data > 0).astype(int)

## merge per Host species and compute prevalence in each host species
heatmap_data = heatmap_data.merge(metadata[[SAMPLE_ID_COL, HOST_COLUMN_NAME]], left_index=True, right_on=SAMPLE_ID_COL)
heatmap_data = heatmap_data.drop(columns=[SAMPLE_ID_COL])
heatmap_data = heatmap_data.groupby(HOST_COLUMN_NAME).sum()

## divide by number of samples per host species
num_samples_per_host = metadata[HOST_COLUMN_NAME].value_counts()  
heatmap_data = heatmap_data.div(num_samples_per_host, axis=0).dropna(axis=0, how='all').T # drop blank rows


# STRIP LABELS AND COLORS
## create strip labels and colors to leaf ID mapping based on most abundant host species per leaf genus
leaf_id_to_host_and_color = {}

for leaf_genus in heatmap_data.index:
    # get host species with highest prevalence for this leaf genus
    host_species = heatmap_data.loc[leaf_genus].idxmax()
    # if two host species have the same prevalence, set label to "Mixed"
    if (heatmap_data.loc[leaf_genus] == heatmap_data.loc[leaf_genus].max()).sum() > 1:
        host_species = "Mixed"
    leaf_id_to_host_and_color[leaf_genus] = (host_species, HOST_COLORS.get(host_species, "lightgrey"))

###########################################################
# CREATE OUTPUT FILES
Path(args.outdir).mkdir(parents=True, exist_ok=True)
## Create the iTOL strip labels and colors files
strip_labels_file = Path(args.outdir) / "leaf_genus_gtdb_to_iTOL_strip_labels.txt"
with open(strip_labels_file, "w") as f:
    f.write("DATASET_COLORSTRIP\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tLeaf_Genus_Host_Species\n")
    f.write("COLOR\t#ff0000\n")
    f.write("COLOR_BRANCHES\t1\n")
    f.write("DATA\n")
    for leaf_id, (host_species, color) in leaf_id_to_host_and_color.items():
        f.write(f"{leaf_id}\t{color}\t{host_species}\n")

## Create the iTOL leaf labels file
leaf_labels_file = Path(args.outdir) / "leaf_genus_gtdb_to_iTOL_labels.txt"
with open(leaf_labels_file, "w") as f:
    f.write("DATASET_LABELS\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tLeaf_Genus_Labels\n")
    f.write("COLOR\t#000000\n")
    f.write("DATA\n")
    for leaf_id, label in leaf_id_to_label.items():
        f.write(f"{leaf_id}\t{label}\n")

## Create the iTOL heatmap file
heatmap_file = Path(args.outdir) / "leaf_genus_gtdb_to_iTOL_heatmap.txt"

## add number of samples per host species to host species names
heatmap_data.columns = [f"{col} (n={num_samples_per_host[col]})" for col in heatmap_data.columns]

with open(heatmap_file, "w") as f:
    f.write("DATASET_HEATMAP\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tLeaf_Genus_Prevalence\n")
    f.write("COLOR\t#ff0000\n")
    f.write("FIELD_LABELS\t" + "\t".join(heatmap_data.columns) + "\n")
    f.write("FIELD_COLORS\t" + "\t".join(["#ff0000"] * len(heatmap_data.columns)) + "\n")
    f.write("DATA\n")
    for leaf_id, row in heatmap_data.iterrows():
        values = "\t".join([f"{v:.3f}" for v in row])
        values = values.replace("0.000", "X") # replace 0.000 with X for itol
        f.write(f"{leaf_id}\t{values}\n")
