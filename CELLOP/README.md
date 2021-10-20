CELLO: Cancer EvoLutionary analysis using LOngitudinal genomic data
================
Jihong TANG
03/10/2021

## Ownership

Wang Lab at HKUST (<http://wang-lab.ust.hk/>)

## Introduction

Cancer EvoLutionary analysis using LOngitudinal genomic data (CELLO) is
a protocal for comprehensive analysis of longitudinal genomic sequencing
data in cancer. It was originally developed by Jiguang Wang \[1,2\], and
was then packed up by Biaobin Jiang, Dong Song and Jihong Tang in MATLAB , R and Python seperately. This code was written in Python.

## Datasets

The input SAVI report (input.savi.txt) consists of a list of genetic
variants from 90 glioblastoma patients.


## Procedure 

```py
python3 code/mutProcess.py -r 20 -a 1 -c 5 \
-g ./example/data/genelist.txt -o ./example/output/ \
./example/data/input.savi.txt
```

```py
python3 code/mutCorrelation.py -c 0.1 -o ./example/figure/ ./example/output/table.mut.gene.txt
```

```py
python3 code/mutFrequency.py -c 5 -g ./example/data/genelist.txt -m ./example/output/table.mut.gene.txt -o ./example/output/ ./example/output/table.savi.txt
```