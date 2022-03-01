#!/bin/bash

source_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/source_data/'
working_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/working'
R_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/Rscripts'

cd $working_dir

## Perform MDS on the merged data

### Extract .genome file from the merged data
plink --bfile merged1 --extract HAS_indepSNP.prune.in --genome --out merged1
#plink --bfile MDS_merge1 --genome --out MDS_merge1

### Perform MDS using PLINK1.9
plink --bfile merged1 --read-genome merged1.genome --cluster --mds-plot 10 --out MDS_merged1
#plink2 --bfile MDS_merge2 --pca 10 --out MDS_merge2_pca
Rscript --no-save $R_dir/MDS_merged.R

### Perform MDS using king
#king -b 1kG_MDS3.bed,HAS_MDS3.bed --mds --cpus 72 --prefix HAS
#Rscripts --no-save $R_dir/MDS_king.R

cd $source_dir
