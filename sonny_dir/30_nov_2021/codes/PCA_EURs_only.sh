#!/bin/bash

## Match the SNPs
awk '{print$2}' HAS_EUR1.bim > HAS_EUR1_SNPs.txt
plink --bfile 1kG_MDS3 --extract HAS_EUR1_SNPs.txt --make-bed --out 1kG_EUR1

## Extract EUR population from the 1kG data
awk '{if ($3 == "CEU" || $3 == "TSI" || $3 == "GBR" || $3 == "FIN") print$1,$2,$3}' pop_1kG.txt > pop_1kG_EUR.txt
awk '{print$1}' pop_1kG_EUR.txt > pop_1kG_EUR_FID.txt
plink --bfile 1kG_EUR1 --keep-fam pop_1kG_EUR_FID.txt --make-bed --out 1kG_EUR2

## Make the combined list of pop labels
awk '{print$1,$2,"HAS"}' HAS_EUR1.fam > pop_HAS_EUR.txt
cat pop_1kG_EUR.txt pop_HAS_EUR.txt | sed -e '1i\IID FID pop' > pop_file_EUR.txt

## Perform PCA on the EUR population
king -b  1kG_EUR2.bed,HAS_EUR1.bed --pca --cpus 72 --prefix PCA_EURs1
