#!/bin/bash

source_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/source_data/'
working_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/working'

cd $working_dir

## 1. Set reference genome
awk '{print$2,$5}' 1kG_MDS2_2.bim > 1kG_ref_list.txt
plink --bfile HAS_MDS2 --reference-allele 1kG_ref_list.txt --make-bed --out HAS_adj

## 2. Resolve strand issues
### 1) Check for possible strand issues
awk '{print$2,$5,$6}' 1kG_MDS2_2.bim > 1kG_MDS2_tmp
awk '{print$2,$5,$6}' HAS_adj.bim > HAS_adj_tmp
sort 1kG_MDS2_tmp HAS_adj_tmp | uniq -u > all_differences.txt

### 2) extract unique entries of the different strands.
awk '{print$1}' all_differences.txt |sort -u > flip_list.txt

### 3) flip the strands in the HAS_adj by taking the 1kG data as the reference.
plink --bfile HAS_adj --flip flip_list.txt --reference-allele 1kG_ref_list.txt --make-bed --out corrected_HAS

### 4) check if the strand issue still persists.
awk '{print$2,$5,$6}' corrected_HAS.bim > corrected_HAS_tmp
sort 1kG_MDS2_tmp corrected_HAS_tmp | uniq -u > uncorresponding_SNPs.txt

### 5) extract SNPs to exclude.
awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exclusion.txt

### 6) exclude the corresopnding SNPs from the HAS and 1kG data
plink --bfile corrected_HAS --exclude SNPs_for_exclusion.txt --make-bed --out HAS_MDS3
plink --bfile 1kG_MDS2_2 --exclude SNPs_for_exclusion.txt --make-bed --out 1kG_MDS3

## 3. Merge HAS and 1kG
plink --bfile HAS_MDS3 --bmerge 1kG_MDS3.bed 1kG_MDS3.bim 1kG_MDS3.fam --make-bed --out merged1


## 4. Convert population codes into superpopulation codes (eg. AFR, AMR, ASN, EUR)
awk '{print$1,$1,$2}' $source_dir/1kG_data/20100804.ALL.panel > pop_1kG.txt
sed 's/JPT/ASN/g' pop_1kG.txt>pop_1kG2.txt # Japanese in Tokyo
sed 's/ASW/AFR/g' pop_1kG2.txt>pop_1kG3.txt # African ancestry in SW USA
sed 's/CEU/EUR/g' pop_1kG3.txt>pop_1kG4.txt  # Utah residents (CEPH) with Northern and Western European ancestry
sed 's/CHB/ASN/g' pop_1kG4.txt>pop_1kG5.txt # Han Chinese
sed 's/CHD/ASN/g' pop_1kG5.txt>pop_1kG6.txt # 
sed 's/YRI/AFR/g' pop_1kG6.txt>pop_1kG7.txt
sed 's/LWK/AFR/g' pop_1kG7.txt>pop_1kG8.txt
sed 's/TSI/EUR/g' pop_1kG8.txt>pop_1kG9.txt #Toscani in Italy
sed 's/MXL/AMR/g' pop_1kG9.txt>pop_1kG10.txt
sed 's/GBR/EUR/g' pop_1kG10.txt>pop_1kG11.txt # British from England and Scotland
sed 's/FIN/EUR/g' pop_1kG11.txt>pop_1kG12.txt # Finnish in Finland
sed 's/CHS/ASN/g' pop_1kG12.txt>pop_1kG13.txt
sed 's/PUR/AMR/g' pop_1kG13.txt>pop_1kG14.txt

## create ethnicity file for the HAS data
awk '{print$1,$2,"HAS"}' HAS_MDS3.fam > pop_HAS.txt
cat pop_1kG14.txt pop_HAS.txt | sed -e '1i\FID IID pop' > pop_file.txt

cd $source_dir
