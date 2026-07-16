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
#           


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

def detect_samples(data_dir, pattern="*"):
        files= glob.glob(os.path.join(data_dir, pattern))
        samples = set()
        
        for f in files:
                base = os.path.basename(f)

                base = re.match(r'^(Pl\d+_[A-Z]\d+)_.*$', base)   # extract sample names
                
                samples.add(base.group(1))
        return sorted(samples)
