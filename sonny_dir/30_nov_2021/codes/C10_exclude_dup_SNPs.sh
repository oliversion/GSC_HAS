#!/bin/bash

source_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/source_data'
working_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/working'
rscripts_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/Rscripts'

cd $working_dir

## check duplicated SNPs in 1kG data
Rscript --no-save $rscripts_dir/check_dup_SNPs.R
## This should produce Routput_dup_SNPs.txt (list of duplicated SNPs from 1kG_MDS1)

## Remove quotation marks from the dup list.
sed 's/\"//g' Routput_dup_SNPs.txt > 1kG_MDS1_dup_SNPs.txt

## Exclude the dup.SNPs from the 1kG data and the HAS data
plink --bfile 1kG_MDS1 --exclude 1kG_MDS1_dup_SNPs.txt --recode --make-bed --out 1kG_MDS2
plink --bfile HAS_MDS1 --exclude 1kG_MDS1_dup_SNPs.txt --make-bed --out HAS_MDS2
### --recode creates .ped and .map file

## Match the builds
awk '{print$2,$4}' HAS_MDS2.map > build_HAS.txt
### build_HAS.txt contains one SNP-id and physical position per line.

plink --bfile 1kG_MDS2 --update-map build_HAS.txt --make-bed --out 1kG_MDS2_2
## 1kG_MDS32 and HAS_MDS2 now have the same build

cd $source_dir
