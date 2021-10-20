#!/usr/bin/env python3

# mutLandscape

import os
import sys
import argparse as ag
import pandas as pd

# import self-defined module
cwd = os.getcwd()
sys.path.append(cwd+'/code/module/')
from landscape import *


def main():
    parser = ag.ArgumentParser(prog='mutProcess', description='mutProcess')
    parser.add_argument('-r', dest='cutoff_refdepth_B',
                        required=True, default=20, type=int, help='')
    parser.add_argument('-a', dest='cutoff_altdepth_B',
                        required=True, default=1, type=int, help='')
    parser.add_argument('-c', dest='cutoff_freq',
                        required=True, default=5, type=int, help='')
    parser.add_argument('-g', dest='genelist', required=True,
                        default='', type=str, help='')
    parser.add_argument('-e', dest='eraselow', required=False,
                        default=True, type=bool, help='')
    parser.add_argument('-o', dest='fp_output', required=False, default='./',
                        type=str, help='Output file path to store the results files. Default: ./')
    parser.add_argument("inputSavi", nargs=1, default='',
                        type=str, help='Input savi table file in txt format')

    # process parameters
    args = parser.parse_args()

    fp_savi = args.inputSavi[0]
    cutoff_refdepth_B = args.cutoff_refdepth_B
    cutoff_altdepth_B = args.cutoff_altdepth_B
    cutoff_freq = args.cutoff_freq
    genelist = args.genelist
    remove_LOW = args.eraselow
    fp_output = args.fp_output

    genelist_selected = pd.read_csv(
        genelist, delimiter='\t', index_col=False, header=0).genelist.to_list()
    df_savi = mR.mutRead(fp_savi, cutoff_refdepth_B,
                         cutoff_altdepth_B, cutoff_freq, fp_output)
    res = mS.mutStats(df_savi, genelist_selected,
                      cutoff_freq, remove_LOW, fp_output)

    if res == 0:
        print("CELLO mutProcess is finished!\n")
    return(res)


if __name__ == '__main__':
    main()
