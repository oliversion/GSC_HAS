#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("snpStats")
library(config)
config = config::get(file = "config_files/config.yml")
library(snpStats)
library(plotly)

#----------Genomic control --------------------------------------------------
data.dir = config$project.dirout.dir = data.dir

# plink --bfile hapmap1 --assoc --adjust --out GC_assoc
plink.GC = read.table(sprintf("%s/GC_assoc.assoc", data.dir), header =TRUE, as.is=T)
median.CHISQ = median(plink.GC$CHISQ)
lambda = median.CHISQ/qchisq(0.5, df=1) # <- 1.09246 PLINK output: 1.08992

hist(plink.GC$CHISQ,main = "Observed chi-squared distribution", xlab = "chi-squared")
