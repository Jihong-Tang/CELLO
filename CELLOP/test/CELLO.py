#!/usr/bin/env python
import sys 
import os
import re
import numpy as np
import pandas as pd
import matplotlib
from scipy.stats import fisher_exact
from scipy.stats import binom_test
# control the working directory
print(os.getcwd())

def main():
    # data import and pandas manipulation
    fp_savi = "./input.savi.txt"
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
    print(cal_mutCorMatrix(df_mutGene, 0.1).head())

    print(mutStats.__doc__)

def mutRead(fp_savi, cutoff_refdepth_B, cutoff_altdepth_B, cutoff_freq):
    """
    Read the savi result table in .txt format and preprocessing the data 
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

    return([df_mut_num, df_mut_gene])


def mutLandscape():
    return(1)


def mutCorrelation(df_mut_gene, cutoff_pValue=0.1):

    # calculate the gene correlation matrix using function cal_mutCorMatrix
    df_mutCorMatrix = cal_mutCorMatrix(df_mut_gene, cutoff_pValue=0.1)

    # plot the Figure 2 as coMutation figure
    return(1)


def mutFrequency(df_savi, genelist_selected, df_mut_gene, cutoff_freq):
    """
    """
    # calculate the mutation frequency matrix
    df_mutFreq = cal_mutFreq(df_savi, genelist_selected, df_mut_gene, cutoff_freq)

    # visualization work
    return(0)

def mutSignature():

    return(0)

def cal_mutCorMatrix(df_mut_gene, cutoff_pValue=0.1):
    """
    Prepare the mutation correlation dataframe for downstream visulization work. 
    Based on the gene based mutation table gained from function mutStats, the 
    correlation between different gene tuples were calculated using Fisher Exact Test.
    Given the cutoff pValue, the gene tuples with significant relationship were filtered
    and were assigned different dot sizes and dot colors based on the pValue and 
    mutation relationship types. The features dot size and dot color are prepared 
    for downstream visulization work and will be illustrated in the docstring of 
    function mutCorrelation().

    Args:
        df_mut_gene (DataFrame): The gene based mutation table in DataFrame format 
            from function mutStats()
        cutoff_pValue (int, default=0.1): The cutoff of the Fisher Exact test pValue
            to filter the relation significance with default value 0.1
    
    Returns:
        DataFrame
    """
    mt_mut_P = df_mut_gene.replace(
        r'[C, P]', 1, regex=True).replace(r'[R, N]', 0, regex=True)
    mt_mut_R = df_mut_gene.replace(
        r'[C, R]', 1, regex=True).replace(r'[P, N]', 0, regex=True)

    # prepare the data of co-mutation analysis
    list_columns = list(df_mut_gene.columns)
    plot_mutCor = [[0] * 4 for _ in range(len(list_columns)**2)]
    idx = 0
    for i in range(len(df_mut_gene.columns)):
        for j in range(len(df_mut_gene.columns)):
            plot_mutCor[idx][0] = list_columns[i]
            plot_mutCor[idx][1] = list_columns[j]

            # initial the matrix for fisher exact test
            mt_FEtest = np.array([[0, 0], [0, 0]])
            if i > j:  # primary tumor section
                mt_FEtest[0][0] = sum((mt_mut_P.iloc[:, i] == 1)
                                      & (mt_mut_P.iloc[:, j] == 1))
                mt_FEtest[0][1] = sum((mt_mut_P.iloc[:, i] == 1)
                                      & (mt_mut_P.iloc[:, j] == 0))
                mt_FEtest[1][0] = sum((mt_mut_P.iloc[:, i] == 0)
                                      & (mt_mut_P.iloc[:, j] == 1))
                mt_FEtest[1][1] = sum((mt_mut_P.iloc[:, i] == 0)
                                      & (mt_mut_P.iloc[:, j] == 0))
                idx_coMutation = ((mt_FEtest[0][0]+1) * (mt_FEtest[1][1]+1)) / \
                    ((mt_FEtest[0][1]+1) * (mt_FEtest[1][0]+1))

                # fisher exact test
                odds, pValue = fisher_exact(mt_FEtest, alternative='two-sided')

                if pValue < cutoff_pValue:
                    # for dot size, if the pvalue is lager than cutoff value, fill in instinct value
                    plot_mutCor[idx][2] = -np.log10(pValue) + 1
                    if idx_coMutation > 1:
                        plot_mutCor[idx][3] = 'D_red' # coMutation 
                    else:
                        plot_mutCor[idx][3] = 'A_blue' # mutual exclusive in primary tumor
                else:
                    plot_mutCor[idx][2] = 1
                    plot_mutCor[idx][3] = 'C_grey'

            elif i < j:  # relapsed tumor section
                mt_FEtest[0][0] = sum((mt_mut_R.iloc[:, i] == 1)
                                      & (mt_mut_R.iloc[:, j] == 1))
                mt_FEtest[0][1] = sum((mt_mut_R.iloc[:, i] == 1)
                                      & (mt_mut_R.iloc[:, j] == 0))
                mt_FEtest[1][0] = sum((mt_mut_R.iloc[:, i] == 0)
                                      & (mt_mut_R.iloc[:, j] == 1))
                mt_FEtest[1][1] = sum((mt_mut_R.iloc[:, i] == 0)
                                      & (mt_mut_R.iloc[:, j] == 0))
                idx_coMutation = ((mt_FEtest[0][0]+1) * (mt_FEtest[1][1]+1)) / \
                    ((mt_FEtest[0][1]+1) * (mt_FEtest[1][0]+1))

                # fisher exact test
                odds, pValue = fisher_exact(mt_FEtest, alternative='two-sided')

                if pValue < cutoff_pValue:
                    # for dot size, if the pvalue is lager than cutoff value, fill in instinct value
                    plot_mutCor[idx][2] = -np.log10(pValue) + 1
                    if idx_coMutation > 1:
                        plot_mutCor[idx][3] = 'E_black' # coMutation
                    else:
                        plot_mutCor[idx][3] = 'B_green' # mutual exclusive in relapsed tumor
                else:
                    plot_mutCor[idx][2] = 1
                    plot_mutCor[idx][3] = 'C_grey'

            else:  # i==j, fill all with 'NA'
                plot_mutCor[idx][2] = np.NaN
                plot_mutCor[idx][3] = np.NaN

            idx += 1  # add the index for plot table

    df_mutCorMatrix = pd.DataFrame(np.array(plot_mutCor), columns=
                                   ['listA', 'listB', 'dotSize', 'dotColor'])
    return(df_mutCorMatrix)


def cal_mutFreq(df_savi, genelist_selected, df_mut_gene, cutoff_freq):
    """
    """
    # remove the synonymous variant labelled by LOW in feature Effect_Impact
    df_savi = df_savi.loc[lambda df: df["Effect_Impact"] != 'LOW', :]
    # pandas Series object ==> list data structure
    case = df_savi["CaseID"].drop_duplicates().tolist()

    # 13 * 3 matrix initial
    mt_mutFreq = [[0] * 3 for _ in range(len(genelist_selected))]
    for i in range(len(genelist_selected)):
        df_savi_gene = df_savi[df_savi.Gene_Name == genelist_selected[i]]
        mt_mutFreq[i][0] = len(df_savi_gene.CaseID[(df_savi_gene.Primary_freq >= cutoff_freq) & (
            df_savi_gene.Recurrent_freq < cutoff_freq)].drop_duplicates())
        mt_mutFreq[i][1] = len(df_savi_gene.CaseID[(df_savi_gene.Primary_freq < cutoff_freq) & (
            df_savi_gene.Recurrent_freq >= cutoff_freq)].drop_duplicates())
        mt_mutFreq[i][2] = sum(df_mut_gene.iloc[:, i] == 'C')

    df_mutFreq = pd.DataFrame(np.array(mt_mutFreq), columns=[
                              'Primary', 'Recurrent ', 'Common'])
    df_mutFreq.index = genelist_selected
    return(df_mutFreq)

def cal_HMDetection(df_savi, cutoff_mut_freq, cutoff_mut_num, cutoff_HM_score):
    return(0)


#def cal_HMFrac()

#def cal_HM_SMratio()
#-------------------------------
if __name__ == "__main__":
    main()
