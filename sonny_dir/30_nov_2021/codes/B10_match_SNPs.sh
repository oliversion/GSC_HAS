#!/bin/bash

source_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/source_data/'
working_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/working'

source_HAS=$working_dir/HAS_init4
source_1kG=$working_dir/1kG_init1

cd $working_dir

## Extract SNP list from the HAS data.
awk '{print$2}' $source_HAS.bim > HAS_SNPs.txt

## Extract the SNPs from the 1kG data corresponding to the list extracted above.
plink --bfile $source_1kG --extract HAS_SNPs.txt --make-bed --out 1kG_MDS1

## Extract SNPs from the 1kG data
awk '{print$2}' 1kG_MDS1.bim > 1kG_MDS1_SNPs.txt

## Extract the SNPs from the HAS data corresponding to the list extracted above.
plink --bfile $source_HAS --extract 1kG_MDS1_SNPs.txt --recode --make-bed --out HAS_MDS1

### Now HAS_MDS1 and 1kG_MDS1 has equal set of SNPs.

cd $source_dir
