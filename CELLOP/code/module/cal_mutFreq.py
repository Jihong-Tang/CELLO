import pandas as pd
import numpy as np

def cal_mutFreq(df_savi, genelist_selected, df_mut_gene, cutoff_freq):
    """
    """
    # remove the synonymous variant labelled by LOW in feature Effect_Impact
    df_savi = df_savi.loc[lambda df: df["Effect_Impact"] != 'LOW', :]
    # pandas Series object ==> list data structure
    # case = df_savi["CaseID"].drop_duplicates().tolist()

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
