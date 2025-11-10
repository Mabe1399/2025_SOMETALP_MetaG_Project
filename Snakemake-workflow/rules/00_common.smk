# COMMON USED FUNCTION SNAKEFILE
# Last Updated: 2025-11-10
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# File containing commonly used functions in the pipelines
# 
# Function included:
#           - detect_samples():
#           - 


##################################################
# IMPORT
##################################################

import os
import glob
import re
import pandas as pd
from snakemake.utils import *

##################################################
# FUNCTIONS
##################################################

def detect_samples(data_dir):
        fastq_files= glob.glob(os.path.join(data_dir,"*.fastq.gz"))
        samples= sorted({re.sub(r"_L\d+_R[12]\.fastq\.gz$","",os.path.basename(f)) for f in fastq_files
        })
        return samples
