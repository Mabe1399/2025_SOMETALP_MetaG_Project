# SOMETALP Metagenomic Pipeline

## Snakemake workflow to get MAGs out of Metagenomice data from multiple hosts

Still under development

### Install 

**SOURCE INSTALL**

Run the commands below:

    git clone https://github.com/Mabe1399/2025_SOMETALP_MetaG_Project.git
    cd 2025_SOMETALP_MetaG_Project

**Once the project is gonna finish the pipeline is gonna cleaned for Reproducibility purposes**

### RUNNING THE CODE
Note: This workflow assumes paired-end multi lane metagenomic data for now

Utilisation is designed to be run locally or in an HPC by following these steps:

1. Fill the profile (local or slurm) file for HPC purposes (designed for slurm scheduler)
   
2. Fill the Config.yaml file for directory specificity and specific tool parameters.

3. Run the snakemake workflow

    ```
    snakemake --profile ~/profile/slurm # or local if run locally
    ```
