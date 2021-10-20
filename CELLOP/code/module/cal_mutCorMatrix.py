import pandas as pd 
import numpy as np 
from scipy.stats import fisher_exact

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