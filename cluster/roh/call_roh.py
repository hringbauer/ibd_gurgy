
'''
Script for Array Job
Runs calling ROH
Runtime: About XX min for 100 individulas
'''

#################################
### Imports  
import socket as socket
import os as os
import sys as sys
import numpy as np
import pandas as pd
from hapsburg.PackagesSupport.hapsburg_run import hapsb_ind 

#################################
### Global Parameters
path_targets = "/mnt/archgen/users/rivollat/Interact/Neo/GRG/20210215_Datafreeze_1240K/genotypes/TF/France_MN_GRG_all_20210215.TF.allChrom"

#################################
### Helper Functions used in Script below

def get_iid(i):
    """Return IID to run.
    Input: Integer index"""
    df2 = pd.read_csv("/mnt/archgen/users/rivollat/Interact/Neo/GRG/20210215_Datafreeze_1240K/genotypes/TF/France_MN_GRG_all_20210215.TF.allChrom.ind", sep="\t", header=None)
    df2.columns = ["iid", "sex", "pop"]

    df21 = pd.read_csv("/mnt/archgen/users/rivollat/Interact/Neo/GRG/20210215_Datafreeze_1240K/genotypes/TF/France_MN_GRG_all_20210215.TF.allChrom.cov.txt", 
                      sep="\t")

    #idx = df21["called SNPs"] > 3e5
    #np.sum(idx)
    df2["n_cov_snps"] = df21["called SNPs"]
    dft = df2[df2["n_cov_snps"]>3e5]
    print(f"Loaded {len(dft)} IIDs with sufficiently many SNPs")
    iids = dft["iid"].values
    return iids[i]


if __name__ == '__main__':  # Only Run if File directly run.
    ### Get the Run Number
    i = int(sys.argv[1]) - 1  # Get the Run number in Python Indexing (qsub)
    print(f"Running Job for Individiual {i}...")
    
    iid = get_iid(i)

    hapsb_ind(iid=iid, chs=range(1, 23), 
              path_targets=path_targets, # The path before the .ind, .snp, .geno
              h5_path1000g='/mnt/archgen/users/hringbauer/data/hapROH.globalRef/chr', 
              meta_path_ref='/mnt/archgen/users/hringbauer/data/hapROH.globalRef/meta_df_all.csv', 
              folder_out="/mnt/archgen/users/hringbauer/git/ibd_gurgy/output/roh/",  # Folder where you want to save the results to 
              processes=1, output=True, delete=False, roh_min_l_final=0.04,
              readcounts=False, logfile=True, combine=True)
    
    print("Job finished! Wuhu!")
    


