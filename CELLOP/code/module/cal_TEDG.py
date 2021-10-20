import pandas as pd
import numpy as np
from scipy.stats import binom_test

def cal_edge(df_mut_gene):
    """
    """

    def which(self):
        try:
            self = list(iter(self))
        except TypeError as e:
            raise Exception("""'which' method can only be applied to iterables.
            {}""".format(str(e)))
        indices = [i for i, x in enumerate(self) if bool(x) == True]
        return(indices)

    df_input = df_mut_gene
    genelist_selected = df_input.columns.tolist()

    # calculating TEDG edge table
    mt_edge = [[0] * len(genelist_selected) for _ in range(len(genelist_selected))]
    list_edge = []

    for i in range(len(genelist_selected)-1):
        start = i+1
        for j in range(start, len(genelist_selected)):
            mt_edge[i][j] = sum((df_input.iloc[:, i] == 'C') &
                                (df_input.iloc[:, j].isin(['P', 'R'])))
            mt_edge[j][i] = sum((df_input.iloc[:, j] == 'C') &
                                (df_input.iloc[:, i].isin(['P', 'R'])))
            # ';'.join() is used to combine the case pattern and is similar to paste function in R
            labelA = ';'.join(df_input.index[which((df_input.iloc[:, i] == 'C') &
                                                (df_input.iloc[:, j].isin(['P', 'R'])))])
            labelB = ';'.join(df_input.index[which((df_input.iloc[:, j] == 'C') &
                                                (df_input.iloc[:, i].isin(['P', 'R'])))])

            if mt_edge[i][j] < mt_edge[j][i]:
                mt_edge[i][j] = 0
                edge = [genelist_selected[j],
                        genelist_selected[i], mt_edge[j][i], labelB]
                list_edge.append(edge)
            elif mt_edge[i][j] > mt_edge[j][i]:
                mt_edge[j][i] = 0
                edge = [genelist_selected[i],
                        genelist_selected[j], mt_edge[i][j], labelA]
                list_edge.append(edge)
            else:
                mt_edge[i][j] = 0
                mt_edge[j][i] = 0
    df_edge = pd.DataFrame(np.array(list_edge), columns=[
                        'geneA', 'geneB', 'weight', 'label'])
    return([df_edge, mt_edge])

def cal_node(df_mut_gene, mt_edge):
    """
    """
    
    df_input = df_mut_gene
    genelist_selected = df_input.columns.tolist()
    mut_freq = np.array([[0] * len(genelist_selected) for _ in range(3)]).T

    for i in range(len(genelist_selected)):
        mut_freq[i][0] = sum(df_input.iloc[:, i] == 'P')
        mut_freq[i][1] = sum(df_input.iloc[:, i] == 'R')
        mut_freq[i][2] = sum(df_input.iloc[:, i] == 'C') * 2

    sample_size = [(mut_freq[i][0] + mut_freq[i][1] + mut_freq[i][2])
                for i in range(len(genelist_selected))]

    tmp_df = pd.DataFrame(np.array(mt_edge))
    ins = [sum(tmp_df.iloc[:, i] > 0) for i in range(len(genelist_selected))]
    outs = [sum(tmp_df.iloc[i, :] > 0) for i in range(len(genelist_selected))]

    pcdf = [1] * len(genelist_selected)
    for i in range(len(pcdf)):
        if ins[i] < outs[i]:
            pcdf[i] = binom_test(x=ins[i], n=ins[i]+outs[i], p=0.5)
        else:
            pcdf[i] = binom_test(x=outs[i], n=ins[i]+outs[i], p=0.5)

    # positive value ==> early; negative value ==> late
    fc = np.log2((np.array(outs)+1) / (np.array(ins)+1))

    tmp_list = [genelist_selected, pcdf, fc, sample_size]
    df_node = pd.DataFrame(np.array(tmp_list).T, columns=[
                        'Gene', 'p_CDF', 'FC', 'Occurence'])
    return(df_node)
