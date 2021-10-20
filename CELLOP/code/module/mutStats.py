import os
import pandas as pd
import numpy as np
def mutStats(df_savi, genelist_selected, cutoff_freq, remove_LOW, fp_output):
    """
    Generate mutation count table and selected gene based mutation type table 
    from the prepocessed savi table and return one list object containing two 
    separate DataFrame objects.
    The mutation count table contains the number of different mutation types (C for 
    Common, P for Primary and R for Recurrent) for each patient with unique caseID. 
    The gene based mutation table contains the mutation type for all the provided
    selected gene for each patient.

    Args:
    df_savi (DataFrame): The savi table in DataFrame format from function mutRead()
    genelist_selected (list): The list of selected gene name to generate the gene
        based mutation type table
    cutoff_freq (int): The cutoff number of mutation frequency to be the mutation 
        detection criterion
    remove_LOW (bool, default=True): The logic tag to determine whether removing 
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

    df_mut_num.to_csv(fp_output+'table.mut.num.txt', index=None, sep='\t')
    df_mut_gene.to_csv(fp_output+'table.mut.gene.txt', index=None, sep='\t')
    return(0)
