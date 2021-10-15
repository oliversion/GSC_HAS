library(config)
config = config::get(file = "config_files/config.yml")
source(sprintf("%s/0_support_func.R", config$project.dir))

# --------------------------------------------------------------------
# ---------------------Data read--------------------------------------
# --------------------------------------------------------------------

project.dir = config$project.dir
data.file = "HA_data_HLA" #"HA_data_1_22"
gwas.fn = lapply(c(bed='bed', bim='bim',fam='fam'), function(n) sprintf("%s/NWE_data/%s.%s", 
                                                                        project.dir, data.file, n))
tic("read.plink")
geno = read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam, na.strings = ("-9"))
toc() # 69.347 sec elapsed

fam = geno$fam

tic("row.summary and col.summary")
snpsum.row <- row.summary(geno$genotypes)
snpsum.col <- col.summary(geno$genotypes)
toc() # 190.029 sec elapsed

# --------------------------------------------------------------------
# 2. ----- Standardized multilocus heterozygosity --------------------
# --------------------------------------------------------------------

tic("smh")
smh = standardized.multilocus.heterozygosity(geno$genotypes, snpsum.row, snpsum.col) 
toc() # 182.509 sec elapsed

# 3. ----------------------- save and load -------------------------------------------------------------
# Write observed and expexted heterozygosity for HA and controls for future use
data.file = sprintf("%s/saved_data/smh_HA_data_1_22_non0MAF.RData", project.dir)
save(fam, snpsum.row, snpsum.col, smh, file = data.file)

tic("data loading")
load(data.file)
toc()

smh_cases    = smh[fam$affected == 2]
smh_controls = smh[fam$affected == 1]

plot.density(fam$affected, smh, project.dir, "SMH", "all_")

# ---------------------------------- Binned standardized.multilocus.heterozygosity ---------------------
smh.binning.mean = sex.binning.stat(smh, mean, geno$fam) 
smh.binning.sd   = sex.binning.stat(smh, sd, geno$fam) 

#----------------------------------- Plot standardized multilocus heterozygosity -----------------------
png(sprintf("%s/plots/SMH/sex_binning_smh.png", project.dir), width = 350, height = 350)
plot.sex.binning.data(smh.binning.mean, smh.binning.sd, "smh")
dev.off() 

plot.stat.cc(smh_cases, smh_controls, project.dir, "SMH", "Standardized multilocus heterozygosity")

#---------------- Mann–Whitney U Test --------------
printf("MAF > 0: %i SNPs\n", dim(snpsum.col)[1])
stat.tests(smh, geno$fam)

#----------------- MAF 0%-1% ---------------------
use = with(snpsum.col, (MAF > 0) & (MAF < 0.01))
printf("MAF [0, 0.01): %i SNPs\n", sum(use))

smh.maf = get.smh(geno$genotypes, use, row.summary, snpsum.col)
stat.tests(smh.maf, geno$fam)

#----------------- MAF 1%-5% ---------------------
use = with(snpsum.col, (MAF >= 0.01) & (MAF < 0.05))
printf("MAF [0.01, 0.05): %i SNPs\n", sum(use))

smh.maf = get.smh(geno$genotypes, use, row.summary, snpsum.col)
stat.tests(smh.maf, geno$fam)

#----------------- MAF 5%-10% ---------------------
use = with(snpsum.col, (MAF >= 0.05) & (MAF < 0.1))
printf("MAF [0.05, 0.1): %i SNPs\n", sum(use))

smh.maf = get.smh(geno$genotypes, use, row.summary, snpsum.col)
stat.tests(smh.maf, geno$fam)

#----------------- MAF 10%-25% ---------------------
use = with(snpsum.col, (MAF >= 0.1) & (MAF < 0.25))
printf("MAF [0.1, 0.25): %i SNPs\n", sum(use))

smh.maf = get.smh(geno$genotypes, use, row.summary, snpsum.col)
stat.tests(smh.maf, geno$fam)

#----------------- MAF 25%-50% ---------------------
use = with(snpsum.col, (MAF >= 0.25) & (MAF <= 0.5))
printf("MAF [0.25, 0.5): %i SNPs\n", sum(use))

smh.maf = get.smh(geno$genotypes, use, row.summary, snpsum.col)
stat.tests(smh.maf, geno$fam)

# --------------------------------------------------------------------
# 2. ---------- Heterozygosity rate ----------------------------------
# --------------------------------------------------------------------

use = with(snpsum.col, Call.rate==1)
printf("SMPs with call.rate =1: %i SNPs\n", sum(use))

hr = get.heterozygosity.rate(geno$genotypes, use)
hr.binning.mean = sex.binning.stat(hr, mean, geno$fam) 
hr.binning.sd   = sex.binning.stat(hr, sd, geno$fam) 
snpsum.row.100 <- row.summary(geno$genotypes[,use])

hr.cases    = hr[geno$fam$affected == 2]
hr.controls = hr[geno$fam$affected == 1]

#----------------------------------- Plot Heterozygosity rate-----------------------
png(sprintf("%s/plots/HR/sex_binning_hr.png", project.dir), width = 350, height = 350)
plot.sex.binning.data(hr.binning.mean, hr.binning.sd, "het rate")
dev.off()

plot.stat.cc(hr.cases, hr.controls, project.dir, "HR", "Heterozygosity rate")

#---------------------------------- TESTS ---------------------------------------------
stat.tests(hr, geno$fam)

#----------------- MAF 0%-1% ---------------------
use = with(snpsum.col, (MAF > 0) & (MAF < 0.01) & (Call.rate==1))
printf("MAF [0, 0.01): %i SNPs\n", sum(use))

hr.maf = get.heterozygosity.rate(geno$genotypes, use)
stat.tests(hr.maf, geno$fam)

#----------------- MAF 1%-5% ---------------------
use = with(snpsum.col, (MAF >= 0.01) & (MAF < 0.05) & (Call.rate==1))
printf("MAF [0.01, 0.05): %i SNPs\n", sum(use))

hr.maf = get.heterozygosity.rate(geno$genotypes, use)
stat.tests(hr.maf, geno$fam)

#----------------- MAF 5%-10% ---------------------
use = with(snpsum.col, (MAF >= 0.05) & (MAF < 0.1) & (Call.rate==1))
printf("MAF [0.05, 0.1): %i SNPs\n", sum(use))

hr.maf = get.heterozygosity.rate(geno$genotypes, use)
stat.tests(hr.maf, geno$fam)

#----------------- MAF 10%-25% ---------------------
use = with(snpsum.col, (MAF >= 0.1) & (MAF < 0.25) & (Call.rate==1))
printf("MAF [0.1, 0.25): %i SNPs\n", sum(use))

hr.maf = get.heterozygosity.rate(geno$genotypes, use)
stat.tests(hr.maf, geno$fam)

#----------------- MAF 25%-50% ---------------------
use = with(snpsum.col, (MAF >= 0.25) & (MAF <= 0.5))
printf("MAF [0.25, 0.5): %i SNPs\n", sum(use))

hr.maf = get.heterozygosity.rate(geno$genotypes, use)
stat.tests(hr.maf, geno$fam)
