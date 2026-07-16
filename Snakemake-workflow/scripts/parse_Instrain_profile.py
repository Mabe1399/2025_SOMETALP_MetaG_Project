#!/usr/bin/env python3
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)
# 
# Usage: python3 parse_Instrain_profile.py --instrain_output instrain_output_dir_list \
# --MAG_representative_info MAG_representative_info_parsed.tsv \
# --stb Scaffold_to_bin_mapping_file.tsv \
# --output_dir output_dir

# output file names:
# Instrain_genome_info_parsed.tsv
# Instrain_gene_info_parsed.tsv
# Instrain_mapping_info_parsed.tsv

# DESCRIPTION
# This script parses the Instrain output files (genome_info.tsv, gene_info.tsv and mapping_info.tsv) for each sample and extracts the relevant information.
# It also adds the sample name to each row and links the scaffolds to bins using the Scaffold_to_bin_mapping_file.tsv file. 
# Finally, it joins the parsed dataframes with the MAG representative info file and adds the corresponding bin name and species name. T
# he parsed dataframes are then saved as tsv files in the output directory.

# IMPORTS
import pandas as pd
import argparse
import os

# PARSE ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument("--instrain_output", required=True, nargs='+', help="List of Instrain output directories")
parser.add_argument("--MAG_representative_info", required=True)
parser.add_argument("--stb", required=True)
parser.add_argument("--output_dir", required=True)
args = parser.parse_args()

# Function to extract genome_info, gene_info, mapping_info and sample_names list from Instrain output directory list
def extract_file_list(instrain_output_dir_list):
    genome_info_files = []
    gene_info_files = []
    mapping_info_files = []

    # extract the sample name from the directory name and use it to identify the corresponding files
    sample_names = [os.path.basename(dir).split('_mapped_')[0] for dir in instrain_output_dir_list]

    # extract the files prefix from dir list and use it to identify the corresponding files
    file_prefixes = [os.path.basename(dir) for dir in instrain_output_dir_list]
    
    # loop through the Instrain output directory list and extract the corresponding files
    for dir in instrain_output_dir_list:
        genome_info_files.append(os.path.join(dir, "output", f"{file_prefixes[instrain_output_dir_list.index(dir)]}_genome_info.tsv"))
        gene_info_files.append(os.path.join(dir, "output", f"{file_prefixes[instrain_output_dir_list.index(dir)]}_gene_info.tsv"))
        mapping_info_files.append(os.path.join(dir, "output", f"{file_prefixes[instrain_output_dir_list.index(dir)]}_mapping_info.tsv"))

    return genome_info_files, gene_info_files, mapping_info_files, sample_names

# Function to parse the genome_info files and add the sample name and extract the relevant information
def parse_genome_info(genome_info_files, sample_names):
    genome_info_df_list = []
    for file in genome_info_files:
        df = pd.read_csv(file, sep='\t')
        df['sample_name'] = sample_names[genome_info_files.index(file)]
        # change header of the first column to "MAG.ID"
        df.rename(columns={df.columns[0]: 'MAG.ID'}, inplace=True)
        # filter per genome breadth >= 0.5 and coverage >= 5
        df = df[(df['breadth'] >= 0.5) & (df['coverage'] >= 5)]
        # calculate actual length of the genome as breadth * length
        df['actual_length'] = df['breadth'] * df['length']
        # calculate normalized read count as read_count / actual_length
        df['normalized_read_count'] = df['filtered_read_pair_count'] / df['actual_length']
        # calculate frequency of Normalized reads mapped as normalized_read_count / sum of normalized_read_count for all genomes in the same sample
        df['frequency_normalized_reads_mapped'] = df['normalized_read_count'] / df['normalized_read_count'].sum()
        genome_info_df_list.append(df)
    return pd.concat(genome_info_df_list, ignore_index=True)

# Function to parse the gene_info files + link scaffolds to bins and add the sample name and extract the relevant information
def parse_gene_info(gene_info_files, sample_names, stb_file):
    gene_info_df_list = []
    stb_df = pd.read_csv(stb_file, sep='\t', header=None, names=['scaffold', 'MAG.ID'])
    for file in gene_info_files:
        df = pd.read_csv(file, sep='\t')
        df['sample_name'] = sample_names[gene_info_files.index(file)]
        df = df.merge(stb_df, left_on='scaffold', right_on='scaffold', how='left')
        gene_info_df_list.append(df)
    return pd.concat(gene_info_df_list, ignore_index=True)

# Function to parse the mapping_info files and extract all Scaffolds rows and add the sample name
def parse_mapping_info(mapping_info_files, sample_names):
    mapping_info_df_list = []
    for file in mapping_info_files:
        df = pd.read_csv(file, sep='\t', skiprows=1)
        df['sample_name'] = sample_names[mapping_info_files.index(file)]
        mapping_info_df_list.append(df[df['scaffold'] == 'all_scaffolds'])
    return pd.concat(mapping_info_df_list, ignore_index=True)

# Function to join any dataframe with the MAG representative info file and add the corresponding bin name and species name
def add_MAG_representative_info(df, MAG_representative_info_file):
    MAG_representative_info_df = pd.read_csv(MAG_representative_info_file, sep='\t')
    df = df.merge(MAG_representative_info_df, left_on='MAG.ID', right_on='MAG ID', how='left')
    return df

# MAIN
if __name__ == "__main__":
    # extract the genome_info, gene_info, mapping_info files and sample names list from the Instrain output directory list
    genome_info_files, gene_info_files, mapping_info_files, sample_names = extract_file_list(args.instrain_output)

    # parse the genome_info files and add the sample name and extract the relevant information
    genome_info_df = parse_genome_info(genome_info_files, sample_names)

    # parse the gene_info files + link scaffolds to bins and add the sample name and extract the relevant information
    gene_info_df = parse_gene_info(gene_info_files, sample_names, args.stb)

    # parse the mapping_info files and extract all Scaffolds rows and add the sample name
    mapping_info_df = parse_mapping_info(mapping_info_files, sample_names)

    # join the genome_info_df with the MAG representative info file and add the corresponding bin name and species name
    genome_info_df = add_MAG_representative_info(genome_info_df, args.MAG_representative_info)

    # join the gene_info_df with the MAG representative info file and add the corresponding bin name and species name
    gene_info_df = add_MAG_representative_info(gene_info_df, args.MAG_representative_info)

    # save the parsed dataframes as tsv files in the output directory
    os.makedirs(args.output_dir, exist_ok=True)
    genome_info_df.to_csv(os.path.join(args.output_dir, "Instrain_genome_info_parsed.tsv"), sep='\t', index=False)
    gene_info_df.to_csv(os.path.join(args.output_dir, "Instrain_gene_info_parsed.tsv"), sep='\t', index=False)
    mapping_info_df.to_csv(os.path.join(args.output_dir, "Instrain_mapping_info_parsed.tsv"), sep='\t', index=False)