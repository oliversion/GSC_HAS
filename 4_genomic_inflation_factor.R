#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("snpStats")
# install.packages("ggplot2")
# install.packages("plotly")

library(snpStats)
library(plotly)

#----------Genomic control --------------------------------------------------
data.dir = paste("/projects/cgstudies/HA_GWAS_2017/Healthy_Aging_GWAS_2017_ABW-P01/Analysis/HAS_Het_Analysis_OLGA/Genomic_contol", sep="")
out.dir = data.dir

# plink --bfile hapmap1 --assoc --adjust --out GC_assoc
plink.GC = read.table(sprintf("%s/GC_assoc.assoc", data.dir), header =TRUE, as.is=T)
median.CHISQ = median(plink.GC$CHISQ)
lambda = median.CHISQ/qchisq(0.5, df=1) # <- 1.09246 PLINK output: 1.08992

hist(plink.GC$CHISQ,main = "Observed chi-squared distribution", xlab = "chi-squared")
