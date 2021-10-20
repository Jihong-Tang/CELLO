#!/usr/bin/env python3

# mutCorrelation

import os
import sys
import argparse as ag
import pandas as pd

# import self-defined module
cwd = os.getcwd()
sys.path.append(cwd+'/code/module/')
import cal_mutCorMatrix as calCM
import plot_mutCorMatrix as pltCM
def main():
    parser = ag.ArgumentParser(prog='mutCorrelation', description='mutCorrelation')
    parser.add_argument('-c', dest='cutoff_pValue', required=False, default=0.1, type=float, help='')
    parser.add_argument('-o', dest='fp_output', required=False, default='./',
                        type=str, help='Output file path to store the result figure. Default: ./')
    parser.add_argument("inputMutGene", nargs=1, default='',
                        type=str, help='Input mutant gene table file in txt format, stored results form mutProcess step')

    # process parameters
    args = parser.parse_args()

    fp_mutgene = args.inputMutGene[0]
    cutoff_pValue = args.cutoff_pValue
    fp_output = args.fp_output

    df_mut_gene = pd.read_csv(fp_mutgene, delimiter='\t')
    df_mutCorMatrix = calCM.cal_mutCorMatrix(df_mut_gene, cutoff_pValue)
    res = pltCM.plot_mutCM(df_mutCorMatrix, fp_output)
    if res == 0:
        print("CELLO mutCorrelation is finished!\n")
    return(res)

if __name__ == '__main__':
    main()
