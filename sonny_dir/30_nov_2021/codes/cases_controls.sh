#!/bin/bash'

awk '{ if ( $6 == "1" ) print $1,$2 }' /projects/cgstudies/HAS_Diversity_OLGA/data/step3_remove_related_and_duplicate_samples/HA_LFS_hapmap_hh_sex_dup_tri_rsID_dup2_strand-flip_tri2_HAonly_maf0_saliva05_related.fam > cases.txt

awk '{ if ($6 == "2") print$1,$2 }' /projects/cgstudies/HAS_Diversity_OLGA/data/step3_remove_related_and_duplicate_samples/HA_LFS_hapmap_hh_sex_dup_tri_rsID_dup2_strand-flip_tri2_HAonly_maf0_saliva05_related.fam > controls.txt

