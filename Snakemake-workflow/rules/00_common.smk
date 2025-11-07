# COMMON USED FUNCTION SNAKEFILE
# Last Updated: 2025-11-07
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# File containing commonly used functions in the pipelines
# 
# Function included:
# get_sample_name_from_Meta(): Extract a list of sample names from the metadata file (.tsv)
#


##################################################
# IMPORT
##################################################

import os
import pandas as pd
from snakemake.utils import *

##################################################
# FUNCTIONS
##################################################

def get_sample_file_path_from_Meta():
    samples = (pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
               .set_index("sample_name", drop=False)
               .sort_index()
               )
    datadir = config.get("raw_data", os.getcwd())
    sample_file_paths = [os.path.join(datadir, name) for name in samples["sample_name"].tolist()]
    
    return sample_file_paths

def get_direction_r1(wildcards) :
    l=get_sample_file_path_from_Meta()

    r_list=[]
    r1_l1=[s for s in l if "L1_R1" in s]
    r1_l2=[s for s in l if "L2_R1" in s]
    r1_l3=[s for s in l if "L3_R1" in s]
    r1_l4=[s for s in l if "L4_R1" in s]

    r_list.extend(r1_l1)
    r_list.extend(r1_l2)
    if len(r1_l3) >0: r_list.extend(r1_l3)
    if len(r1_l4) >0: r_list.extend(r1_l4)
    
    return(r_list)

def get_direction_r2(wildcards) :
    l=get_sample_file_path_from_Meta()

    r_list=[]
    r2_l1=[s for s in l if "L1_R2" in s]
    r2_l2=[s for s in l if "L2_R2" in s]
    r2_l3=[s for s in l if "L3_R2" in s]
    r2_l4=[s for s in l if "L4_R2" in s]

    r_list.extend(r2_l1)
    r_list.extend(r2_l2)
    if len(r2_l3) >0: r_list.extend(r2_l3)
    if len(r2_l4) >0: r_list.extend(r2_l4)
    
    return(r_list)

