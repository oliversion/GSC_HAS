## Load .bim data from 1kG
setwd("/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/working")
mds_data = read.table("1kG_MDS1.bim", sep="\t")

## Check duplicated rows
dat = mds_data[duplicated(mds_data$V2)==1,]
dups = data.frame(V1= as.character(dat$V2))

## Export the duplicated SNPs
write.table(dups, "Routput_dup_SNPs.txt", col.names = F, row.names = F)

