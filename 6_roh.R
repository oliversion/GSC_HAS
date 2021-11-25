# import libraries====================================================================
library(config)
library(snpStats)
library(tictoc)
library(ggplot2)
config = config::get(file = "roh_config.yml")
source("0_support_func.R")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("snpStats")

#====================================================================
#      READ DATA
# ====================================================================
project.dir = config$project.dir

# plink results for runs of homozygosity
plink.roh = read.table(sprintf("%s/%s", config$project.dir, config$plink.hom.indiv), header =TRUE, as.is=T)
plink.roh_autosomal = read.table(sprintf("%s/%s", config$project.dir, config$plink.hom.indiv.autosomal), header =TRUE, as.is=T)
plink.roh_x = read.table(sprintf("%s/%s", config$project.dir, config$plink.hom.indiv.X), header =TRUE, as.is=T)

# plink results for calculating three inbreeding coefficients for each sample
# Fhat1	Variance-standardized relationship minus 1
# Fhat2	Excess homozygosity-based inbreeding estimate (same as PLINK --het)
# Fhat3	Estimate based on correlation between uniting gametes
plink.ibc = read.table(sprintf("%s/%s", config$project.dir, config$plink.ibc), header =TRUE, as.is=T)

# PCs after preprocessing made by Rawnak
plink.pcs = read.table(sprintf("%s/%s", config$project.dir, config$plink.pcs), header =TRUE, as.is=T)
plink.pcs = plink.pcs[, 1:12]

# raw autosomal data
data.file = config$plink.autosomal
gwas.fn = lapply(c(bed='bed', bim='bim',fam='fam'), function(n) sprintf("%s/%s.%s",project.dir, data.file, n))
tic("data readin...")
geno = read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam, na.strings = ("-9"))
toc() # 58.062 sec elapsed

geno.fam = geno$fam[,c('member', 'sex', 'affected')]
# ('1' = male, '2' = female, '0' = unknown)
geno.fam$sex[geno.fam$sex==1]='male'
geno.fam$sex[geno.fam$sex==2]='female'
geno.fam$IID=geno.fam$member
f = formula("PHE~FROH + sex+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")

#====================================================================
#      TEST RUNS OF HOMOZYGOSITY FOR AUTOSOMAL
# ====================================================================
df.merged = base::merge(x = plink.roh_autosomal, y = plink.pcs, by =c('IID'), all.x = TRUE)
df.merged = base::merge(x = df.merged, y = geno.fam[,c('IID', 'sex')], by =c('IID'), all.x = TRUE)
df.merged$PHE[df.merged$PHE==1]=0
df.merged$PHE[df.merged$PHE==2]=1
df.merged$PHE = as.factor(df.merged$PHE)
df.merged$sex = as.factor(df.merged$sex)
df.merged$FROH = df.merged$KB/(3e+09)

# summary statistics
with(df.merged, tapply(KB, PHE, mean))
# 0        1 
# 10085.66 11215.95 
with(df.merged, tapply(KB, PHE, median))
# 0        1 
# 8989.960 8825.245 

# use logistic regression to test
fit = glm(f, family = binomial, data = df.merged)
s = summary(fit)
print(s)

p1 = ggplot(df.merged, aes(x = sex, y = KB, color=PHE)) + 
  geom_boxplot() + 
  ggtitle("FROH for autosomal chr")

#====================================================================
#      TEST RUNS OF HOMOZYGOSITY FOR X CHR
# ====================================================================
df.merged = base::merge(x = plink.roh_x, y = plink.pcs, by =c('IID'), all.x = TRUE)
df.merged = base::merge(x = df.merged, y = geno.fam[,c('IID', 'sex')], by =c('IID'), all.x = TRUE)
df.merged$PHE[df.merged$PHE==1]=0
df.merged$PHE[df.merged$PHE==2]=1
df.merged$PHE = as.factor(df.merged$PHE)
df.merged$sex = as.factor(df.merged$sex)
df.merged$FROH = df.merged$KB/(3e+09)
fit = glm(f, family = binomial, data = df.merged)
s = summary(fit)
print(s)
p1 = ggplot(df.merged, aes(x = sex, y = FROH, color=PHE)) + 
  geom_boxplot() + 
  ggtitle("FROH for X chr")

#====================================================================
#      TEST RUNS OF HOMOZYGOSITY FOR ALL
# ====================================================================
df.merged = base::merge(x = plink.roh, y = plink.pcs, by =c('IID'), all.x = TRUE)
df.merged = base::merge(x = df.merged, y = geno.fam[,c('IID', 'sex')], by =c('IID'), all.x = TRUE)
df.merged$PHE[df.merged$PHE==1]=0
df.merged$PHE[df.merged$PHE==2]=1
df.merged$PHE = as.factor(df.merged$PHE)
df.merged$sex = as.factor(df.merged$sex)
df.merged$FROH = df.merged$KB/(3e+09)
fit = glm(f, family = binomial, data = df.merged)
s = summary(fit)
print(s)

#====================================================================
#      TEST F3 FOR AUTOSOMAL
# ====================================================================
f.f3 = formula("PHE~Fhat3 + sex+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")

df.merged = base::merge(x = plink.ibc, y = plink.pcs, all.x = TRUE)
df.merged = base::merge(x = df.merged, y = geno.fam[,c('IID', 'sex', 'affected')], by =c('IID'), all.x = TRUE)
df.merged$PHE = df.merged$affected
df.merged$PHE[df.merged$PHE==1]=0
df.merged$PHE[df.merged$PHE==2]=1
df.merged$PHE = as.factor(df.merged$PHE)
df.merged$sex = as.factor(df.merged$sex)

# summary statistics
with(df.merged, tapply(Fhat3, PHE, mean))
# 0               1 
# 0.0005662947 0.0011434841 
with(df.merged, tapply(Fhat3, PHE, median))
# 0               1 
# -4.86686e-05  4.24250e-04  

fit = glm(f.f3, family = binomial, data = df.merged)
s = summary(fit)
print(s)

p1 = ggplot(df.merged, aes(x = sex, y = Fhat3, color=PHE)) + 
  geom_boxplot() + 
  ggtitle("F3 for autosomal chr")

p1 = ggplot(df.merged, aes(x = PHE, y = Fhat3)) +
  geom_violin(aes(color = PHE), trim = FALSE, position = position_dodge(0.9), fill = NA) +
  geom_boxplot(width = 0.2) +
  facet_grid( ~ sex) +
   ylab("F3") +
  xlab("Phenotype")+
  theme(legend.position = "none")+
  ggtitle("Results for F3, p-value = 0.0124")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        aspect.ratio=1)
