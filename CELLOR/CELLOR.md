CELLO: Cancer EvoLutionary analysis using LOngitudinal genomic data
================
SONG, Dong
1/15/2020

## Ownership

Wang Lab at HKUST (<http://wang-lab.ust.hk/>)

## Introduction

Cancer EvoLutionary analysis using LOngitudinal genomic data (CELLO) is
a protocal for comprehensive analysis of longitudinal genomic sequencing
data in cancer. It was originally developed by Jiguang Wang \[1,2\], and
was then packed up by Biaobin Jiang and Dong Song in MATLAB and R
seperately. This code was written in R.

## Datasets

The input SAVI report (input.savi.txt) consists of a list of genetic
variants from 90 glioblastoma patients.

## Loading required R packages

``` r
library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
library(ggtern)
```

    ## Registered S3 methods overwritten by 'ggtern':
    ##   method           from   
    ##   +.gg             ggplot2
    ##   grid.draw.ggplot ggplot2
    ##   plot.ggplot      ggplot2
    ##   print.ggplot     ggplot2

    ## --
    ## Remember to cite, run citation(package = 'ggtern') for further info.
    ## --

    ## 
    ## Attaching package: 'ggtern'

    ## The following objects are masked from 'package:gridExtra':
    ## 
    ##     arrangeGrob, grid.arrange

    ## The following objects are masked from 'package:ggplot2':
    ## 
    ##     %+%, aes, annotate, calc_element, ggplot, ggplot_build,
    ##     ggplot_gtable, ggplotGrob, ggsave, layer_data, theme, theme_bw,
    ##     theme_classic, theme_dark, theme_gray, theme_light, theme_linedraw,
    ##     theme_minimal, theme_void

``` r
library(ggalt)
```

    ## Registered S3 methods overwritten by 'ggalt':
    ##   method                  from   
    ##   grid.draw.absoluteGrob  ggplot2
    ##   grobHeight.absoluteGrob ggplot2
    ##   grobWidth.absoluteGrob  ggplot2
    ##   grobX.absoluteGrob      ggplot2
    ##   grobY.absoluteGrob      ggplot2

``` r
library("igraph")
```

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

## Loading CELLO functions

……

## CELLO Pipeline

``` r
savi.table<-somaticfilter("../../input.savi.txt",20,1,5)
head(savi.table[,1:12])
```

    ##     chr       pos ref alt                              Effect Effect_Impact
    ## 1 chr10   5248362   T   C splice_donor_variant+intron_variant          HIGH
    ## 2 chr10  27703023   A   T                    missense_variant      MODERATE
    ## 3 chr10 122263428   G   A                    missense_variant      MODERATE
    ## 4 chr11    209598   C   T                  synonymous_variant           LOW
    ## 5 chr11  20676292  TG   T                  frameshift_variant          HIGH
    ## 6 chr12  21422680   C   T                  synonymous_variant           LOW
    ##   Functional_Class Codon_Change Amino_Acid_Change Amino_Acid_length Gene_Name
    ## 1                -            -                 -               323    AKR1C4
    ## 2         MISSENSE      Tcc/Acc              S53T               767    PTCHD3
    ## 3         MISSENSE      cGc/cAc              R52H               271  PPAPDC1A
    ## 4           SILENT      tcC/tcT              S108               537     RIC8A
    ## 5                -         tgg/              W758               797    SLC6A5
    ## 6           SILENT      ccG/ccA              P605               670   SLCO1A2
    ##   Sgt1_max_frequency
    ## 1                 35
    ## 2                 32
    ## 3                 42
    ## 4                 12
    ## 5                 21
    ## 6                 35

``` r
knownDriverGene <- c('LTBP4','PTPN11','NF1','RB1','PDGFRA','PIK3CG','PIK3R1','PIK3CA','PTEN','EGFR','IDH1','ATRX','TP53')

mutNum.table <- mutStatistics(savi.table,5)

head(mutNum.table)
```

    ##   Patients Primary Common Recurrent
    ## 1     R001       4     56         2
    ## 2     R002      23     46        19
    ## 3     R003      49      1         1
    ## 4     R004       2     45         4
    ## 5     R005      23     50        16
    ## 6     R006       5      6       165

``` r
mutGenes.table <- mutGenes(savi.table, knownDriverGene,5,remove_LOW = TRUE)
head(mutGenes.table)
```

    ##      LTBP4 PTPN11 NF1 RB1 PDGFRA PIK3CG PIK3R1 PIK3CA PTEN EGFR IDH1 ATRX TP53
    ## R001 "N"   "N"    "N" "N" "N"    "N"    "N"    "C"    "N"  "N"  "N"  "N"  "C" 
    ## R002 "N"   "N"    "C" "N" "N"    "C"    "C"    "N"    "N"  "N"  "N"  "N"  "N" 
    ## R003 "N"   "N"    "N" "N" "N"    "N"    "N"    "N"    "N"  "N"  "N"  "N"  "N" 
    ## R004 "N"   "N"    "C" "N" "N"    "N"    "N"    "C"    "N"  "N"  "N"  "N"  "N" 
    ## R005 "N"   "N"    "N" "N" "N"    "N"    "N"    "C"    "N"  "P"  "N"  "N"  "N" 
    ## R006 "R"   "N"    "R" "N" "N"    "N"    "N"    "R"    "N"  "N"  "N"  "N"  "R"

``` r
mutLandscape(mutNum.table,mutGenes.table)
```

<img src="Output/img/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

``` r
coMutation(mutGenes.table)
```

    ## Warning: Removed 13 rows containing missing values (geom_point).

<img src="Output/img/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

``` r
freq.table <- freqMutation(savi.table, knownDriverGene,mutGenes.table,5)
```

<img src="Output/img/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

``` r
HM.table <- hyperMutation(savi.table,15,350,1.2)
```

<img src="Output/img/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

    ## P-value between Primary and NonHM Recurrence: 0.229260285986123
    ## P-value between Primary and HM Recurrence: 2.59649100398582e-06
    ## P-value between NonHM Recurrence and HM Recurrence: 4.38254782225301e-05

``` r
Cluster.table <- evoCluster(mutNum.table)
```

    ## Warning: Solution to limits produces range outside of [0,1] for some scales

<img src="Output/img/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

``` r
switch.table <- mutSwitch(savi.table,knownDriverGene,5,20)
```

<img src="Output/img/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

``` r
selGene <-c('LTBP4','IDH1','ATRX','TP53','NF1','MSH6','PIK3CG','PIK3R1','PIK3CA','PTEN','EGFR')
allMutGenes.table <- mutGenes(savi.table, selGene,5,remove_LOW = TRUE)
TEDGedge.table <- getTEDG(allMutGenes.table)
```

<img src="Output/img/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

``` r
TEDGedge.table
```

    ##    geneA    geneB    weight label                               
    ## 1  "IDH1"   "LTBP4"  "3"    "R027;R049;R055"                    
    ## 2  "ATRX"   "LTBP4"  "1"    "R049"                              
    ## 3  "TP53"   "LTBP4"  "3"    "R027;R049;R055"                    
    ## 4  "PTEN"   "LTBP4"  "1"    "R039"                              
    ## 5  "EGFR"   "LTBP4"  "2"    "R022;R039"                         
    ## 6  "IDH1"   "TP53"   "1"    "R054"                              
    ## 7  "IDH1"   "NF1"    "5"    "R043;R049;R051;R053;R055"          
    ## 8  "IDH1"   "MSH6"   "3"    "R027;R044;R051"                    
    ## 9  "IDH1"   "PIK3R1" "1"    "R051"                              
    ## 10 "IDH1"   "PIK3CA" "2"    "R048;R049"                         
    ## 11 "IDH1"   "PTEN"   "2"    "R046;R051"                         
    ## 12 "IDH1"   "EGFR"   "1"    "R027"                              
    ## 13 "ATRX"   "TP53"   "1"    "R054"                              
    ## 14 "ATRX"   "NF1"    "5"    "R043;R049;R051;R053;R077"          
    ## 15 "ATRX"   "MSH6"   "2"    "R044;R051"                         
    ## 16 "ATRX"   "PIK3R1" "1"    "R051"                              
    ## 17 "ATRX"   "PIK3CA" "2"    "R048;R049"                         
    ## 18 "EGFR"   "ATRX"   "1"    "R039"                              
    ## 19 "TP53"   "NF1"    "7"    "R034;R043;R049;R051;R053;R055;R077"
    ## 20 "TP53"   "MSH6"   "3"    "R027;R044;R051"                    
    ## 21 "TP53"   "PIK3CG" "1"    "R042"                              
    ## 22 "TP53"   "PIK3CA" "2"    "R048;R049"                         
    ## 23 "TP53"   "PTEN"   "5"    "R034;R038;R046;R051;R061"          
    ## 24 "EGFR"   "TP53"   "2"    "R022;R039"                         
    ## 25 "NF1"    "PIK3CG" "1"    "R042"                              
    ## 26 "PIK3R1" "NF1"    "1"    "R100"                              
    ## 27 "PTEN"   "MSH6"   "1"    "R039"                              
    ## 28 "EGFR"   "MSH6"   "2"    "R022;R039"                         
    ## 29 "PTEN"   "PIK3CG" "1"    "R042"                              
    ## 30 "PIK3R1" "PTEN"   "1"    "R024"

## Reference

\[1\] Wang, J., Cazzato, E., Ladewig, E., Frattini, V., Rosenbloom, D.
I., Zairis, S., … & Lee, J. K. (2016). Clonal evolution of glioblastoma
under therapy. Nature Genetics, 48(7), 768-776.

\[2\] Wang, J., Khiabanian, H., Rossi, D., Fabbri, G., Gattei, V.,
Forconi, F., … & Pasqualucci, L. (2014). Tumor evolutionary directed
graphs and the history of chronic lymphocytic leukemia. Elife, 3,
e02869.
