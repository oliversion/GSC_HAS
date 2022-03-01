#!/bin/bash

## Compute GC of the original data (HAS_MDS3)
plink --bfile HAS_MDS3 --maf 0.01 --hwe 0.01 --exclude /projects/cgstudies/HAS_Diversity_OLGA/data/annotation/inversion.txt --range --indep-pairwise 50 5 0.5 --out HAS_indepSNP3

plink --bfile HAS_MDS3 --extract HAS_indepSNP3.prune.in --make-bed --out HAS_MDS3_pruned

plink --bfile HAS_MDS3_pruned --assoc --adjust --out HAS_MDS3_assoc



## Apply filters on the EUR pop data

### Extract SNPs that satisty: MAF 0.01, HWE 0.01, exclude HLA regeion, LDpruning 50 5 0.5
#plink --bfile HAS_EUR1 --maf 0.01 --hwe 0.01 --exclude /projects/cgstudies/HAS_Diversity_OLGA/data/annotation/inversion.txt --range --indep-pairwise 50 5 0.5 --out HAS_indepSNP2

### Make pruned plink file
#plink --bfile HAS_EUR1 --extract HAS_indepSNP2.prune.in --make-bed --out HAS_EUR1_pruned

### Compute the Genomic inflation factor score: value closer to 1, the better.
#plink --bfile HAS_EUR1_pruned --assoc --adjust --out HAS_EUR1_assoc
