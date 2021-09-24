#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import matplotlib
# control the working directory
print(os.getcwd())

def main():
    # data import and pandas manipulation
    fp_savi = "../../input.savi.txt"
    cutoff_refdepth_B = 20
    cutoff_altdepth_B = 1
    cutoff_freq = 5
    df_savi = mutRead(fp_savi, cutoff_refdepth_B, cutoff_altdepth_B, cutoff_freq)
    print(df_savi.head())

def mutRead(fp_savi, cutoff_refdepth_B, cutoff_altdepth_B, cutoff_freq):
    """Read the savi result table in .txt format and preprocessing the data 
    based on several cutoffs, returning a DataFrame object storing all the 
    filtered information.

    Args:
    fp_savi (str): The filepath of the savi result file to read in
    cutoff_refdepth_B (int): The cutoff number of feature refdepth_Blood to
        preprocess the data
    cutoff_altdepth_B (int): The cutoff number of feature altdepth_Blood to 
        preprocess the data
    cutoff_freq (int): The cutoff number of feature Sgt1_max_frequency to 
        preprocess the data

    Returns:
    DataFrame
    """
    raw_data = pd.read_csv(fp_savi, delimiter="\t")
    df_savi = raw_data.loc[lambda df:(df["refdepth_Blood"] >= cutoff_refdepth_B) &
                                 (df["altdepth_Blood"] <= cutoff_altdepth_B) &
                                 (df["Sgt1_max_frequency"] >= cutoff_freq), :]
    df_savi = df_savi.fillna(0)
    return(df_savi)

def mutStats(df_savi ):
    return(df_savi)

#-------------------------------
if __name__ == "__main__":
    main()
