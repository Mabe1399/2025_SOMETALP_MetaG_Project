#!/usr/bin/env python3
# Author: Matias Becker Burgos (Matias.B)
# Script to subset metagenomic data based on number of reads or total base pairs
# Usage: python Subset_data.py -i input.fastq.gz -o output.fastq.gz -n 1000000
# 

# Imports
import gzip
import random
import argparse

# Function to subset FASTQ files
def subset_fastq(fastq1, fastq2, output1, output2, num_reads, seed=None):
    if seed is not None:
        random.seed(seed)
    
    # count total reads
    print("Counting total reads in input files...")
    with gzip.open(fastq1, 'rt') as f1:
        total_reads = sum(1 for _ in f1) // 4
    print(f"Total reads in input: {total_reads}")

    if num_reads > total_reads:
        raise ValueError("Requested number of reads exceeds total reads in input files.")
    
    # Reservoir sampling
    keep_indices = set(random.sample(range(total_reads), num_reads))
    print(f"Selected {num_reads} reads to keep.")

    # Write subsetted reads to output files
    kept_reads =0
    with gzip.open(fastq1, 'rt') as f1, gzip.open(fastq2, 'rt') as f2, \
            gzip.open(output1, 'wt') as out1, gzip.open(output2, 'wt') as out2:
            
            for i in range(total_reads):
                read1 = [next(f1) for _ in range(4)]
                read2 = [next(f2) for _ in range(4)]
                
                if i in keep_indices:
                    out1.writelines(read1)
                    out2.writelines(read2)
                    kept_reads += 1
    print(f"Wrote {kept_reads} reads to output files.")

# Main function to parse arguments and call subset function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subset paired FASTQ files based on number of reads.")
    parser.add_argument("-i1", "--input1", required=True, help="Input FASTQ R1 (.gz).")
    parser.add_argument("-i2", "--input2", required=True, help="Input FASTQ R2 (.gz).")
    parser.add_argument("-o1", "--output1", required=True, help="Output FASTQ R1 (.gz).")
    parser.add_argument("-o2", "--output2", required=True, help="Output FASTQ R2 (.gz).")
    parser.add_argument("-n", "--num_reads", type=int, required=True, help="Number of reads to subset.")
    parser.add_argument("-s", "--seed", type=int, default=None, help="Random seed")
    
    args = parser.parse_args()
    
    subset_fastq(args.input1, args.input2, args.output1, args.output2, args.num_reads, args.seed)