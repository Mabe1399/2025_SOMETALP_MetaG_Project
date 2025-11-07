# QUALITY CONTROL AND PREPROCESSING SNAKEFILE
# Last Updated: 2025-11-07
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# QC and Preprocessing rules usable in many different pipelines
# The workflow start with the assumption that all the data files are in the same directory 
# and the sequencing occured in many lanes.

# Steps included:
# 01: Concatenation (When multiple lane sequencing)
# 02: Before QC (FastQC + MultiQC)
# 03: Trimming (Trimmomatic)
# 04: After QC (FastQC + MultiQC)
# 05: Before Taxonomic profiling (Kraken 2)
# 06: Host filtering (BWA)
# 07: After Taxonomic profiling (Kraken 2)

##################################################
# CONFIGURATION
##################################################
configfile: "../../config/config.yaml"

# If `working_dir` is set in the config file, tell Snakemake to use it.
# Use a safe default (the current working directory) when the key is missing.
workdir: config.get("working_dir", "../../.")

##################################################
# SET UP VARIABLES
##################################################

import os
import pandas as pd
from snakemake.utils import *


##################################################
# RULE ALL
##################################################

#rule all:


##################################################
# RULES
##################################################

############ 01 Concat Lanes #####################



