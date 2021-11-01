# import libraries====================================================================
library(config)
library(snpStats)
library(tictoc)
config = config::get(file = "roh_config.yml")
source("0_support_func.R")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("snpStats")

#====================================================================
#      READ DATA
# ====================================================================
project.dir = config$project.dir

plink.roh = read.table(sprintf("%s/%s", config$project.dir, config$plink.hom.indiv), header =TRUE, as.is=T)
plink.roh_autosomal = read.table(sprintf("%s/%s", config$project.dir, config$plink.hom.indiv.autosomal), header =TRUE, as.is=T)
plink.roh_x = read.table(sprintf("%s/%s", config$project.dir, config$plink.hom.indiv.X), header =TRUE, as.is=T)

plink.pcs = read.table(sprintf("%s/%s", config$project.dir, config$plink.pcs), header =TRUE, as.is=T)
plink.pcs = plink.pcs[, 1:12]

# raw autosomal data
data.file = config$plink.autosomal
gwas.fn = lapply(c(bed='bed', bim='bim',fam='fam'), function(n) sprintf("%s/%s.%s",project.dir, data.file, n))
tic("data readin...")
geno = read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam, na.strings = ("-9"))
toc() # 58.062 sec elapsed

geno.fam = geno$fam[,c('member', 'sex')]
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
#      PLOTS
# ====================================================================
p1 = ggplot(df.merged, aes(x = sex, y = FROH, color=PHE)) + 
  geom_boxplot() + 
  ggtitle("FROH")
