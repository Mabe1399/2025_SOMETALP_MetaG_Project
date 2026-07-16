#!/usr/bin/env python3
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)
# Script to combine multiple Kraken2 reports into a single report
# Usage: python Combine_kraken_reports.py -i report1.txt report2.txt ... -o combined_report.txt

# output format: taxid name rank sample1_pct sample1_taxa_reads sample1_clade_reads sample2_pct sample2_taxa_reads sample2_clade_reads ...

# Imports
import argparse
import pandas as pd 
from pathlib import Path


# Function to combine Kraken2 reports
def parse_kraken_reports(file_path):
    # Extract sample name from file path
    sample_name = Path(file_path).stem.replace('_kraken2', '')

    # Read the Kraken2 report into a DataFrame
    df = pd.read_csv(
        file_path, 
        sep='\t', 
        header=None, 
        names=[
        'pct', 'reads_clade', 'reads_taxon', 'rank', 'taxid', 'name'],
        dtype={'taxid': str}
        )
    
    # Select relevant columns and rename them
    df = df[['taxid', 'name', 'rank', 'pct', 'reads_taxon', 'reads_clade']]
    df = df.rename(columns={
        "pct": f'{sample_name}_pct',
        "reads_taxon":  f'{sample_name}_taxa_reads',
        "reads_clade": f'{sample_name}_clade_reads'
        })
    
    return df

def main():
    # Argument parser
    parser = argparse.ArgumentParser(description='Combine multiple Kraken2 reports into a single report.')
    parser.add_argument('-i', '--input', nargs='+', required=True, help='Input Kraken2 report files.')
    parser.add_argument('-o', '--output', required=True, help='Output combined report file.')
    args = parser.parse_args()

    # Initialize an empty DataFrame for the combined report
    combined_df = pd.DataFrame()

    # Process each input file
    for file_path in args.input:
        df = parse_kraken_reports(file_path)
        
        if combined_df.empty:
            combined_df = df
        else:
            combined_df = pd.merge(combined_df, df, on=['taxid', 'name', 'rank'], how='outer')

    # Sort the combined DataFrame by taxid
    combined_df = combined_df.sort_values(by='taxid', key=lambda x: x.astype(int))

    # Save the combined report to the output file
    combined_df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()