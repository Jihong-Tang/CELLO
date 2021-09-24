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

    genelist_selected = ['LTBP4', 'PTPN11', 'NF1', 'RB1', 'PDGFRA',
                     'PIK3CG', 'PIK3R1', 'PIK3CA', 'PTEN', 'EGFR', 'IDH1', 'ATRX', 'TP53']
    cutoff_freq = 5
    remove_LOW = True
    [df_mutNum, df_mutGene] = mutStats(df_savi, genelist_selected, 5, remove_LOW)
    print(df_mutNum.head())
    print(df_mutGene.head())

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


def mutStats(df_savi, genelist_selected, cutoff_freq, remove_LOW=True):
    """Generate mutation count table and selected gene based mutation type table 
    from the prepocessed savi table and return one list object containing two 
    separate DataFrame objects.
    The mutation count table contains the number of different mutation types (C for 
    Common, P for Primary and R for Recurrent) for each patient with unique caseID. 
    The gene based mutation table contains the mutation type for all the provided
    selected gene for each patient.

    Args:
    df_savi (DataFrame): The savi table in DataFrame from function mutRead()
    genelist_selected (list): The list of selected gene name to generate the gene
        based mutation type table
    cutoff_freq (int): The cutoff number of mutation frequency to be the mutation 
        detection criterion
    remove_LOW (bool, default = True): The logic tag to determine whether removing 
        the synonymous variant labelled by LOW in feature Effect_Impact is necessary

    Returns:
    list(DataFrame, DataFrame) 
    """
    # remove the synonymous variant labelled by LOW in feature Effect_Impact
    if(remove_LOW):
        df_savi = df_savi.loc[lambda df: df["Effect_Impact"] != 'LOW', :]
    
    # pandas Series object ==> list data structure
    case = df_savi["CaseID"].drop_duplicates().tolist()

    # initial the storing list to count different mutations
    mut_P = [0] * len(case)
    mut_R = [0] * len(case)
    mut_C = [0] * len(case)

    # use lambda function to count mutations in different type
    for i in range(len(case)):
        mut_C[i] = (lambda df: (df["CaseID"] == case[i])
                    & (df["Primary_freq"] >= cutoff_freq)
                    & (df["Recurrent_freq"] >= cutoff_freq)
                    )(df_savi).sum()
        mut_P[i] = (lambda df: (df["CaseID"] == case[i])
                    & (df["Primary_freq"] >= cutoff_freq)
                    & (df["Recurrent_freq"] < cutoff_freq)
                    )(df_savi).sum()
        mut_R[i] = (lambda df: (df["CaseID"] == case[i])
                    & (df["Primary_freq"] < cutoff_freq)
                    & (df["Recurrent_freq"] >= cutoff_freq)
                    )(df_savi).sum()

    df_mut_num = pd.DataFrame({'Patient': case,
                                'Primary': mut_P,
                                'Common': mut_C,
                                'Recurrent': mut_R})
    
    # create 2D list in python using following format to make sure to create several
    # separate lists in the second dimension
    list_mut_gene = [['N'] * len(genelist_selected) for _ in range(len(case))]
    df_savi_gene = df_savi[df_savi.Gene_Name.isin(genelist_selected)]

    geneSel = genelist_selected
    for i in range(len(case)):
        for j in range(len(geneSel)):
            temp_P = df_savi_gene.Primary_freq[(df_savi_gene.CaseID == case[i]) & (
                df_savi_gene.Gene_Name == geneSel[j])]
            temp_R = df_savi_gene.Recurrent_freq[(df_savi_gene.CaseID == case[i]) & (
                df_savi_gene.Gene_Name == geneSel[j])]
            if any(temp_P >= cutoff_freq) & any(temp_R >= cutoff_freq):
                list_mut_gene[i][j] = 'C'
            elif any(temp_P >= cutoff_freq):
                list_mut_gene[i][j] = 'P'
            elif any(temp_R >= cutoff_freq):
                list_mut_gene[i][j] = 'R'
            else:
                list_mut_gene[i][j] = 'N'
    
    df_mut_gene = pd.DataFrame(np.array(list_mut_gene), columns=geneSel)
    df_mut_gene.index = case

    return([df_mut_num, df_mut_gene])
#-------------------------------
if __name__ == "__main__":
    main()
