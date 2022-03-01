#!/bin/bash

target_data='/projects/cgstudies/HAS_Diversity_OLGA/data/step3_remove_related_and_duplicate_samples/HA_LFS_hapmap_hh_sex_dup_tri_rsID_dup2_strand-flip_tri2_HAonly_maf0_saliva05_related'

wwrking_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/working'
source_dir='/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/source_data'

##### QC on 1000 Genomes data #####

cd $working_dir

plink --bfile $target_data --make-bed --out HAS_init0

## Extract autosomal SNPs from the HAS data
#awk '{if ($1 >= 1 && $1 <=22) print$2}' HAS_init.bim > HAS_autosomes.txt
#plink --bfile HAS_init --extract HAS_autosomes.txt --make-bed --out HAS_init1
plink --bfile HAS_init0 --chr 1-22 --make-bed --out HAS_init1

## Remove SNPs with missing genotype rate > 3%
plink --bfile HAS_init1 --geno 0.03 --make-bed --out HAS_init2

## Remove SNPs with HWE p-values <1x10-6
plink --bfile HAS_init2 --hwe 1e-6 --make-bed --out HAS_init3

## Remove SNPs with MAF < 0.01 (1%) 
plink --bfile HAS_init3 --maf 0.01 --make-bed --out HAS_init4

## LD pruning
plink --bfile HAS_init4 --indep-pairwise 50 5 0.5 --out HAS_indepSNP


##### QC on 1000 Genomes data #####

plink --bfile $source_dir/1kG_data/ALL.2of4intersection.20100804.genotypes_no_missing_IDs --allow-no-sex --make-bed --out $working_dir/1kG_init

## Extract SNPs from autosomes only
#awk '{ if($1>=1 && $1<=22) print$2}' 1kG_init.bim > 1kG_autosomes.txt
#plink --bfile 1kG_init --extract 1kG_autosomes.txt --make-bed --out 1kG_init1
plink --bfile 1kG_init --chr 1-22 --make-bed --out 1kG_init1

cd $source_dir
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
