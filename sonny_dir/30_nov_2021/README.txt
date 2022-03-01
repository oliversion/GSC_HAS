#### Ethnicity checks using PLINK and KING ####

## Here, we perform two dimension reduction methods: PCA and MDS

## To run MDS, run A10 to E10 .sh 
## To run PCA, run A10 to D10 then E11 .sh

## Please read below for the descriptions and codes for each file.


# ------------------------------------------------------------------------------------

# A10_qc_HAS_1kG.sh

# We weren't sure which data in the 'data' folder we should use as our super-senior data. 
# So we picked one that looked fairly reasonable.
# The thresholds for MAF, HWE p-values, and the missing genotype rate referred to the values on /projects/cgstudies/HAS_Diversity_OLGA/data/plink_qc.sh

target_data=target_data='/projects/cgstudies/HAS_Diversity_OLGA/data/step3_remove_related_and_duplicate_samples/HA_LFS_hapmap_hh_sex_dup_tri_rsID_dup2_strand-flip_tri2_HAonly_maf0_saliva05_related' # HAS data dir

source_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/source_data' # 1kG data dir

working_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/working' ## Sonnys working dir

rscripts_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/Rscripts' ## Sonnys rscript dir



## Copy target_data and save it as HAS_init0
plink --bfile $target_data --make-bed --out HAS_init0

## Extract autosomal SNPs from the HAS data
plink --bfile HAS_init0 --chr 1-22 --make-bed --out HAS_init1

## Remove SNPs with missing genotype rate > 3%
plink --bfile HAS_init1 --geno 0.03 --make-bed --out HAS_init2

## Remove SNPs with HWE p-values <1x10-6
plink --bfile HAS_init2 --hwe 1e-6 --make-bed --out HAS_init3

## Remove SNPs with MAF < 0.005 (0.05%) 
plink --bfile HAS_init3 --maf 0.005 --make-bed --out HAS_init4

## LD pruning
plink --bfile HAS_init4 --indep-pairwise 50 5 0.5 --out HAS_indepSNP

## QC on 1000 Genomes data #####

plink --bfile $source_dir/1kG_data/ALL.2of4intersection.20100804.genotypes_no_missing_IDs --allow-no-sex --make-bed --out $working_dir/1kG_init  # 1kG data downloaded on 

## Extract SNPs from autosomes only
plink --bfile 1kG_init --chr 1-22 --make-bed --out 1kG_init1

# END
# -----------------------------------------------------------------------------------

# B10_match_SNPs.sh

source_HAS=$working_dir/HAS_init4
source_1kG=$working_dir/1kG_init1


## Extract SNP list from the HAS data.
awk '{print$2}' $source_HAS.bim > HAS_SNPs.txt

## Extract the SNPs from the 1kG data corresponding to the list extracted above.
plink --bfile $source_1kG --extract HAS_SNPs.txt --make-bed --out 1kG_MDS1

## Extract SNPs from the 1kG data
awk '{print$2}' 1kG_MDS1.bim > 1kG_MDS1_SNPs.txt

## Extract the SNPs from the HAS data corresponding to the list extracted above.
plink --bfile $source_HAS --extract 1kG_MDS1_SNPs.txt --recode --make-bed --out HAS_MDS1

### Now HAS_MDS1 and 1kG_MDS1 has equal set of SNPs.

# END
# ------------------------------------------------------------------------------------

# C10_exclude_dup_SNPs.sh


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

# END
# -------------------------------------------------------------------------------------

# D10_merge_HAS_1kG.sh

# Here we merge HAS and 1kG, and then create ethnicity data based on the merged data.

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
sed 's/TSI/EUR/g' pop_1kG8.txt>pop_1kG9.txt # Toscani in Italy
sed 's/MXL/AMR/g' pop_1kG9.txt>pop_1kG10.txt
sed 's/GBR/EUR/g' pop_1kG10.txt>pop_1kG11.txt # British from England and Scotland
sed 's/FIN/EUR/g' pop_1kG11.txt>pop_1kG12.txt # Finnish in Finland
sed 's/CHS/ASN/g' pop_1kG12.txt>pop_1kG13.txt
sed 's/PUR/AMR/g' pop_1kG13.txt>pop_1kG14.txt

## create ethnicity file for the HAS data
awk '{print$1,$2,"HAS"}' HAS_MDS3.fam > pop_HAS.txt
cat pop_1kG14.txt pop_HAS.txt | sed -e '1i\FID IID pop' > pop_file.txt


# END 
# ---------------------------------------------------------------------------------------------

# This is where we perform the dimensionality reduction methods.

# E10_perform_MDS.sh ---> perform MDS using PLINK1.9 and KING.

## Extract .genome file from the merged data
plink --bfile merged1 --extract HAS_indepSNP.prune.in --genome --out merged1

## Perform MDS using PLINK1.9
plink --bfile merged1 --read-genome merged1.genome --cluster --mds-plot 10 --out MDS_merged1
Rscript --no-save $R_dir/MDS_merged.R

## Perform MDS using king
king -b 1kG_MDS3.bed,HAS_MDS3.bed --mds --cpus 72 --prefix HAS
Rscripts --no-save $rscripts_dir/MDS_king.R


# E11_perform_PCA.sh ---> perform PCA using PLINK2.0 and KING.

## Extract .genome file from the merged data
plink --bfile merged1 --extract HAS_indepSNP.prune.in --genome --out merged1

## Perform PCA using PLINK 2.0
plink2 --bfile merged1 --pca 10 --out PCA_merged1  ## We did not know how to run PCA on PLINK1.9. We used PLINK2.0 instead.
Rscript --no-save $R_dir/PCA_merged.R  ## creates PCA plots (PC1 vs. PC2, PC2 vs. PC3, PC1 vs. PC3) in .pdf
 
### Perform PCA using KING
king -b 1kG_MDS3.bed,HAS_MDS3.bed --pca --cpus 72 --prefix HAS 
Rscripts --no-save $R_dir/PCA_king.R  ## creates PCA plots (PC1 vs. PC2, PC2 vs. PC3, PC1 vs. PC3) in .pdf

# END
# ------------------------------------------------------------------------------------------------







