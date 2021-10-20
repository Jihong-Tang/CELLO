#!/usr/bin/env python3

# mutSwitch

import os
import sys
import argparse as ag
import pandas as pd

# import self-defined module
cwd = os.getcwd()
sys.path.append(cwd+'/code/module/')
import cal_mutSwitch as calmS
import plot_mutSwitch as pltmS

def main():
    parser = ag.ArgumentParser(
        prog='mutCorrelation', description='mutCorrelation')
    parser.add_argument('-c', dest='cutoff_freq',
                        required=True, default=5, type=int, help='')
    parser.add_argument('-g', dest='genelist', required=True,
                        default='', type=str, help='')
    parser.add_argument('-m', dest='fp_mutgene', required=True, default='',
                        type=str, help='')
    parser.add_argument('-o', dest='fp_output', required=False, default='./',
                        type=str, help='Output file path to store the result figure. Default: ./')
    parser.add_argument("inputSavi", nargs=1, default='',
                        type=str, help='Input mutant gene table file in txt format, stored results form mutProcess step')

    # process parameters
    args = parser.parse_args()

    fp_savi = args.inputSavi[0]
    cutoff_freq = args.cutoff_freq
    genelist = args.genelist
    fp_mutgene = args.fp_mutgene
    fp_output = args.fp_output

    df_mut_gene = pd.read_csv(fp_mutgene, delimiter='\t')
    df_savi = pd.read_csv(fp_savi, delimiter='\t')
    genelist_selected = pd.read_csv(
        genelist, delimiter='\t', index_col=False, header=0).genelist.to_list()
    df_mutFreq = calFq.cal_mutFreq(
        df_savi, genelist_selected, df_mut_gene, cutoff_freq)
    print(df_mutFreq.head())
    #if res == 0:
    #    print("CELLO mutCorrelation is finished!\n")
    #return(res)


if __name__ == '__main__':
    main()
