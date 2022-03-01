#!/bin/bash

source_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/source_data/'
working_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/working'
R_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/Rscripts'
cd $working_dir

## create ethnicity file for the HAS data
awk '{print$1,$2,"HAS"}' HAS_MDS3.fam > pop_HAS.txt
cat pop_1kG14.txt pop_HAS.txt | sed -e '1i\FID IID pop' > pop_file.txt

## Perform MDS on the merged data

### Extract .genome file from the merged data
plink --bfile merged1 --extract HAS_indepSNP.prune.in --genome --out merged1
#plink --bfile MDS_merge1 --genome --out MDS_merge1

### Perform PCA using PLINK 2.0
#plink2 --bfile merged1 --pca 10 --out PCA_merged1
Rscript --no-save $R_dir/PCA_merged.R
 
### Perform PCA using KING
king -b 1kG_MDS3.bed,HAS_MDS3.bed --pca --cpus 72 --prefix HAS 
Rscripts --no-save $R_dir/PCA_king.R

cd $source_dir
