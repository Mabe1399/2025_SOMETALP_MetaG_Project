#!/usr/bin/env python3
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)
# 
# Usage: python3 parse_dRep_genomeInfo.py --CheckM_summary CheckM_summary_file.tsv --output output_file.tsv

# Output file: dRep_genome_info_parsed.tsv

# IMPORTS
import pandas as pd
import argparse

# PARSE ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument("--CheckM_summary", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

# LOAD CheckM SUMMARY
checkm_summary = pd.read_csv(args.CheckM_summary, sep="\t")

# EXTRACT GENOME INFO FROM CheckM SUMMARY
## Extract genome ID, completeness and contamination

genome_info = checkm_summary[["Bin Id", "Completeness", "Contamination"]]
genome_info.columns = ["genome", "completeness", "contamination"]

# Add .fasta extension to genome ID
genome_info["genome"] = genome_info["genome"].apply(lambda x: x + ".fa")

# SAVE OUTPUT
genome_info.to_csv(args.output, sep=",", index=False)