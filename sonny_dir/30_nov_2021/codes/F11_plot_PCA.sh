#!/bin/bash

source_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/source_data'
working_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/working'
rscripts_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/Rscripts'

cd $working_dir

## convert population codes into superpopulation codes (eg. AFR, AMR, ASN, EUR)
awk '{print$1,$1,$2}' $source_dir/1kG_data/20100804.ALL.panel > pop_1kG.txt
sed 's/JPT/ASN/g' pop_1kG.txt>pop_1kG2.txt
sed 's/ASW/AFR/g' pop_1kG2.txt>pop_1kG3.txt
sed 's/CEU/EUR/g' pop_1kG3.txt>pop_1kG4.txt 
sed 's/CHB/ASN/g' pop_1kG4.txt>pop_1kG5.txt
sed 's/CHD/ASN/g' pop_1kG5.txt>pop_1kG6.txt
sed 's/YRI/AFR/g' pop_1kG6.txt>pop_1kG7.txt
sed 's/LWK/AFR/g' pop_1kG7.txt>pop_1kG8.txt
sed 's/TSI/EUR/g' pop_1kG8.txt>pop_1kG9.txt
sed 's/MXL/AMR/g' pop_1kG9.txt>pop_1kG10.txt
sed 's/GBR/EUR/g' pop_1kG10.txt>pop_1kG11.txt
sed 's/FIN/EUR/g' pop_1kG11.txt>pop_1kG12.txt
sed 's/CHS/ASN/g' pop_1kG12.txt>pop_1kG13.txt
sed 's/PUR/AMR/g' pop_1kG13.txt>pop_1kG14.txt

## create ethnicity file for the HAS data
awk '{print$1,$2,"HAS"}' HAS_MDS3.fam > pop_HAS.txt
cat pop_1kG14.txt pop_HAS.txt | sed -e '1i\FID IID pop' > pop_file.txt

## Plot PC1 vs. PC2
#Rscript --no-save $rscripts_dir/PCA_merged.R

#cp PCA_result.pdf ~

cd $source_dir
