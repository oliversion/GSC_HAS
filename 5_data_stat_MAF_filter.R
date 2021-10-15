library(snpStats)
library(tictoc)
library(config)
config = config::get(file = "config_files/config.yml")

# ---------------------Data read------------------------------------------------------------------------
project.dir = config$project.dir
source(sprintf("%s/support_func.R", project.dir))

data.file = "HA_data_1_22"
gwas.fn = lapply(c(bed='bed', bim='bim',fam='fam'), function(n) sprintf("%s/NWE_data/%s.%s", 
                                                                        project.dir, data.file, n))
geno = read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam, na.strings = ("-9"))

snpsum.row <- row.summary(geno$genotypes)
snpsum.col <- col.summary(geno$genotypes)

# Call rate
png(sprintf("%s/plots/HA_filters_as_stat.png", project.dir), width = 350, height = 350)
boxplot(snpsum.row$Call.rate, snpsum.col$Call.rate,
        main = "Filtered autosomal SNP statistic (cr>99%, pHWE<5e-5)",
        # main = "NWE data autosomal statistics",
        at = c(1,2),
        names = c("Sample call rate", "SNP call rate")
)
dev.off()

# MAF distribution
png(sprintf("%s/plots/HA_data_1_22_MAF.png", project.dir), width = 350, height = 350)
hist(snpsum.col$MAF,main = "MAF distribution", xlab = "MAF")
dev.off()
#plink.imiss= read.table(sprintf("%s/plink_output/autosome_cr.imiss", project.dir), header =TRUE, as.is=T) 

#--------------------------------------------- X CHR ------------------------------------
is.female = geno$fam$sex == 2
geno.female = geno$genotypes[is.female]
female.snpsum.row <- row.summary(geno.female)
female.snpsum.col <- col.summary(geno.female)
png(sprintf("%s/plots/NWE_X_stat.png", project.dir), width = 350, height = 350)
boxplot(female.snpsum.row$Call.rate, female.snpsum.col$Call.rate,
        main = "NWE data X chr statistics",
        at = c(1,2),
        names = c("Sample call rate", "SNP call rate")
)
dev.off()

#--------------------------------------------- MAF ------------------------------------
# Combined MAF
# plink --bfile HA_data_1_22 --freq --noweb --out plink_output/MAF_1_22
plink.MAF = read.table(sprintf("%s/plink_output/MAF_1_22.frq", project.dir), header =TRUE, as.is=T)

# MAF for case/controls
#plink --bfile HA_data_1_22 --freq case-control --noweb --out plink_output/MAF_cc
plink.MAF.cc = read.table(sprintf("%s/plink_output/MAF_cc.frq.cc", project.dir), header =TRUE, as.is=T)
maf.cases = plink.MAF.cc$MAF_A
maf.controls = plink.MAF.cc$MAF_U

# MAF for case/controls female/male 
# plink --bfile HA_data_1_22 --freq case-control --filter-males --out plink_output/MAF_cc_male
plink.MAF.cc.female = read.table(sprintf("%s/plink_output/MAF_cc_female.frq.cc", project.dir), header =TRUE, as.is=T)
plink.MAF.cc.male   = read.table(sprintf("%s/plink_output/MAF_cc_male.frq.cc", project.dir), header =TRUE, as.is=T)

snp.use = get.non0.maf(plink.MAF.cc, plink.MAF.cc.female, plink.MAF.cc.male)
snp.to.exclude = plink.MAF[!snp.use,]
plink.maf.not.0 = plink.MAF[snp.use,] # 2 681 695 SNPs

write.table(snp.to.exclude$SNP, file = sprintf("%s/NWE_data/zero_MAF.txt", project.dir), sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote=FALSE)

#---------------------------------------MAF X CHR ------------------------------------
# plink --bfile NWE_data/HA_data_X --freq case-control --out plink_output/MAF_X_cc
plink.MAF.cc = read.table(sprintf("%s/plink_output/MAF_X_cc.frq.cc", project.dir), header =TRUE, as.is=T)
maf.cases = plink.MAF.cc$MAF_A
maf.controls = plink.MAF.cc$MAF_U
snp.to.exclude = (plink.MAF.cc$MAF_U == 0)|(plink.MAF.cc$MAF_A == 0)
write.table(plink.MAF.cc[snp.to.exclude, "SNP"], file = sprintf("%s/NWE_data/zero_MAF_X.txt", project.dir), sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote=FALSE)

#---------------------------------------MAF X CHR for equal GP ancestry ------------------------------------
# plink --bfile NWE_data/HA_data_X_gp_equal --freq case-control --out plink_output/MAF_X_gp_equal
plink.MAF.cc = read.table(sprintf("%s/plink_output/MAF_X_gp_equal.frq.cc", project.dir), header =TRUE, as.is=T)
maf.cases = plink.MAF.cc$MAF_A
maf.controls = plink.MAF.cc$MAF_U
snp.to.exclude = (plink.MAF.cc$MAF_U == 0)|(plink.MAF.cc$MAF_A == 0)
write.table(plink.MAF.cc[snp.to.exclude, "SNP"], file = sprintf("%s/NWE_data/zero_MAF_X_gp_equal.txt", project.dir), sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote=FALSE)

#---------------------------------------MAF 1-22 for equal GP ancestry ------------------------------------
# plink --bfile NWE_data/HA_data_X --freq case-control --out plink_output/MAF_X_cc
plink.MAF.cc = read.table(sprintf("%s/plink_output/MAF_1_22_gp_equal.frq.cc", project.dir), header =TRUE, as.is=T)
maf.cases = plink.MAF.cc$MAF_A
maf.controls = plink.MAF.cc$MAF_U
snp.to.exclude = (plink.MAF.cc$MAF_U == 0)|(plink.MAF.cc$MAF_A == 0)
write.table(plink.MAF.cc[snp.to.exclude, "SNP"], file = sprintf("%s/NWE_data/zero_MAF_1_22_gp_equal.txt", project.dir), sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote=FALSE)

#--------------------------------------plot-------------------------------------------
png(sprintf("%s/plots/maf_cc.png", project.dir), width = 350, height = 350)
plot(maf.cases, maf.controls, main = "MAF HA vs controls",xlab = "HA", ylab = "controls", pch='.')
dev.off() 

png(sprintf("%s/plots/HA_data_1_22_MAF_plink.png", project.dir), width = 350, height = 350)
hist(plink.MAF$MAF,main = "MAF distribution plink", xlab = "MAF")
dev.off()

