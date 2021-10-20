import pandas as pd
import numpy as np
import re

def cal_mutSwitch(df_savi, genelist_selected, cutoff_freq_low, cutoff_freq_high):
    """
    """
    # define new which function from R to pandas, code coming from
    # https://alex.miller.im/posts/python-pandas-which-function-indices-similar-to-R/

    def which(self):
        try:
            self = list(iter(self))
        except TypeError as e:
            raise Exception("""'which' method can only be applied to iterables.
            {}""".format(str(e)))
        indices = [i for i, x in enumerate(self) if bool(x) == True]
        return(indices)

    df_savi = df_savi.loc[lambda df: df["Effect_Impact"] != 'LOW', :]
    case = df_savi.CaseID.drop_duplicates().tolist()

    caseGene = []
    mutType = []
    x_content = []
    y_content = []
    color_factor = []

    for i in range(len(case)):  # for each patient
        tmp_case = df_savi[df_savi.CaseID == case[i]]

        for j in range(len(genelist_selected)):
            nameCom = case[i] + '-' + genelist_selected[j]
            tmp_switch = tmp_case[tmp_case.Gene_Name == genelist_selected[j]]
            if any((tmp_switch.Primary_freq <= cutoff_freq_low) & (tmp_switch.Recurrent_freq >= cutoff_freq_high)):
                if any((tmp_switch.Primary_freq >= cutoff_freq_high) & (tmp_switch.Recurrent_freq <= cutoff_freq_low)):
                    switch_list = which(((tmp_switch.Primary_freq <= cutoff_freq_low) & (tmp_switch.Recurrent_freq >= cutoff_freq_high)) | (
                        (tmp_switch.Primary_freq >= cutoff_freq_high) & (tmp_switch.Recurrent_freq <= cutoff_freq_low)))
                    for k in switch_list:
                        caseGene += [nameCom, nameCom]
                        AAC = re.sub(
                            '^ *| *$', '', tmp_switch.Amino_Acid_Change.tolist()[k])
                        if AAC in mutType:
                            tmp = AAC + 'NA'
                            mutType += [tmp, tmp]
                        else:
                            mutType += [AAC, AAC]
                        x_content += ['Primary', 'Recurrence']
                        y_content += [tmp_switch.Primary_freq.tolist()[k], tmp_switch.Recurrent_freq.tolist()[k]]
                        if tmp_switch.Primary_freq.tolist()[k] > tmp_switch.Recurrent_freq.tolist()[k]:
                            color_factor += ['P', 'P']
                        else:
                            color_factor += ['R', 'R']

    tmp_list = [caseGene, mutType, x_content, y_content, color_factor]
    df_switch = pd.DataFrame(np.array(tmp_list).T, columns=['CaseGene', 'MutType', 'x_content', 'y_content', 'color_factor'])
    return(df_switch)
    #if len(df_switch) == 0:
    #    print('No switch event detected in the given gene list!\n')
