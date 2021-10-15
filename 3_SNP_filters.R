# BiocManager::install("VariantAnnotation")
# BiocManager::install("ensemblVEP") to parse VEP CSQ field

library(tictoc)
library(VariantAnnotation)
library(ensemblVEP)
library(snpStats)
library(config)
config = config::get(file = "config_files/config.yml")

# https://bioconductor.org/packages/release/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf

project.dir = config$project.dir

#Original data in plink format
data.file = "NWE_final_cohort" # HA_data_HLA: filtered for 6:25.5M-33.5M
gwas.fn = lapply(c(bed='bed', bim='bim',fam='fam'), function(n) sprintf("%s/NWE_data/%s.%s",project.dir, data.file, n))
geno = read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam, na.strings = ("-9"))

# ====================================================================
#            POSITION MAP
# ====================================================================
# save SNP with posiion for further analysis
snp.map = geno$map[, c("chromosome", "snp.name", "position")]
colnames(snp.map) = c("CHR", "SNP", "position")

write.table(snp.map, file = sprintf("%s/NWE_data/NWE_final_cohort_map.snp", project.dir), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)

# ====================================================================
#            SNP RANGES
# ====================================================================

snp.map$start = snp.map$position-3
snp.map$end = snp.map$position+2 
#using bedtools anf fasta file: starts with 0, does not include the last position

snp.map$CHR[snp.map$CHR==23] = "X"
snp.map$CHR[snp.map$CHR==24] = "Y"
snp.map$CHR[snp.map$CHR==25] = "X"
snp.map$CHR[snp.map$CHR==26] = "MT"

options("scipen"=100, "digits"=4) #a penalty for scientific display
write.table(snp.map[,c("CHR", "start", "end")], file = sprintf("%s/NWE_data/NWE_final_cohort_map_range.bed", project.dir), sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote=FALSE)

# ====================================================================
#            FAMILY BACKGROUND ANALYSIS
# ====================================================================
#read query
HA.query = read.table(sprintf("%s/NWE_data/HAQuery_191126_processed_csv.csv", project.dir), 
                      header =TRUE, sep=",", stringsAsFactors=FALSE)
HA.db.samples = HA.query[,c("Sample.Name..", "MatGF_Ancestry", "MatGM_Ancestry",
                            "PatGF_Ancestry", "PatGM_Ancestry")]
colnames(HA.db.samples) = c("member", "MatGF_Ancestry", "MatGM_Ancestry",
                            "PatGF_Ancestry", "PatGM_Ancestry")
#setdiff(samples.in.study$member,res$member)
#mismatch in sample id "139_ROB"  "294_MCI"  "C496_GIL" "C509_SAN" "C596_JOH"

samples.in.study = geno$fam[,c("pedigree", "member", "affected")]
samples.in.study = samples.in.study[with(samples.in.study, order(member)),]
res = base::merge(x = samples.in.study, y = HA.db.samples, by="member", x.all = TRUE)

is.match = res$MatGF_Ancestry==res$MatGM_Ancestry &  res$MatGF_Ancestry==res$PatGF_Ancestry &
            res$MatGF_Ancestry==res$PatGM_Ancestry
matches.gp = res[is.match,]
not.matches.gp = res[!is.match,]

write.table(not.matches.gp[not.matches.gp$affected==1,], file = "not_matched_controls", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)
write.table(not.matches.gp[not.matches.gp$affected==2,], file = "not_matched_HA", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)

write.table(matches.gp[, c("pedigree", "member")], file = "matched_samples", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote=FALSE)
# ====================================================================
#            HLA filters
# ====================================================================
# HLA Region: MHC6:28477797-33448354
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25#/def_region-MHC
snp.filtered = geno$map[geno$map$chromosome == 6 & geno$map$position>=28477797 & geno$map$position<=33448354,
                        c("chromosome", "snp.name")]
colnames(snp.filtered) = c("CHR", "SNP")
write.table(snp.filtered[c("CHR", "SNP")], file = sprintf("%s/NWE_data/annotation/HLA.snp", project.dir), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)

# HLA extended: Region: MHC6:25.5M-33.5M
snp.filtered = geno$map[geno$map$chromosome == 6 & geno$map$position>=25500000 & geno$map$position<=33500000,
                         c("chromosome", "snp.name")]
colnames(snp.filtered) = c("CHR", "SNP")
write.table(snp.filtered[c("CHR", "SNP")], file = sprintf("%s/NWE_data/annotation/HLA_ext.snp", project.dir), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)

snp.list = geno$map$snp.name

# ====================================================================
#            LCT filters
# ====================================================================
#https://ghr.nlm.nih.gov/gene/LCT#location
# 135,787,850 to 135,837,195
snp.filtered = geno$map[geno$map$chromosome == 2 & geno$map$position>=135100000 & geno$map$position<=136800000,
                        c("chromosome", "snp.name")]
colnames(snp.filtered) = c("CHR", "SNP")
write.table(snp.filtered[c("CHR", "SNP")], file = sprintf("%s/NWE_data/annotation/LCT_ext.snp", project.dir), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)

# ====================================================================
#            dbSNP annotation in vcf format
# ====================================================================
vcf.file = "NWE_data/annotation/dbSNP/dbSNP_annotation_filtering.vcf"

#check for info columns  
hdr = scanVcfHeader(vcf.file)
info.header = info(hdr)
info.fields = rownames(info.header)

# info field that contains dbSNP Database codes: R3 R5 SYN NSN NSM NSF U3 U5 INT ASS DSS REF
info.list = c("R3", "R5", "SYN", "NSN", "NSM", "NSF", "U3", "U5", "INT", "ASS", "DSS", "REF")
param = ScanVcfParam(info=info.list)

tic("vcf file loading")
vcf = readVcf(vcf.file, param=param)
toc()

vcf.rows = rowRanges(vcf)
df.CHR = data.frame(seqnames(vcf.rows))
colnames(df.CHR) = "CHR"
df.SNP = data.frame(ranges(vcf.rows))
df.SNP = df.SNP[c("names", "start")]
colnames(df.SNP) = c("SNP", "position")

# vcf.fixed = fixed(vcf)
vcf.info = info(vcf)

dbSNP.df = cbind(df.CHR,df.SNP, vcf.info)
head(dbSNP.df)

dbSNP.snp.list = dbSNP.df$SNP

# ====================================================================
#           SNP without annotation
# ====================================================================
snp.without.annotation = setdiff(snp.list,dbSNP.snp.list) #4002 SNPs

write.table(snp.without.annotation, file = sprintf("%s/NWE_data/annotation/snp_without_annotation.snp", project.dir), sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote=FALSE)

# ====================================================================
#     Variant Effect Predictor results (for snp.without.annotation)
# ====================================================================
vcf.vep.file = "NWE_data/annotation/VEP/VEP_for_non_in_dbSNP_vcf.vcf"
vcf.vep = readVcf(vcf.vep.file) # CSQ field
csq.rows =  parseCSQToGRanges(vcf.vep, VCFRowID=rownames(vcf.vep))
df.CHR = data.frame(seqnames(csq.rows))
colnames(df.CHR) = "CHR"
df.SNP = data.frame(ranges(csq.rows))
df.SNP = df.SNP[c("names", "start")]
colnames(df.SNP) = c("SNP", "position")

vep.Consequence = csq.rows$Consequence
# head(VEP.df)


# ====================================================================
#            Saving filters
# ====================================================================

#1. Intronic SNP: FxnCode = 6; INFO Field = INT
snp.filtered = dbSNP.df[dbSNP.df$INT == TRUE,]
write.table(snp.filtered[c("CHR", "SNP")], file = sprintf("%s/NWE_data/annotation/intronic.snp", project.dir), sep = "\t", 
           row.names = FALSE, col.names = TRUE, quote=FALSE)

#2. UTR SNP: FxnCode = 53,55; INFO Field = U3, U5
snp.filtered = dbSNP.df[dbSNP.df$U3 == TRUE | dbSNP.df$U5 == TRUE,]
write.table(snp.filtered[c("CHR", "SNP")], file = sprintf("%s/NWE_data/annotation/UTR.snp", project.dir), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)

#3. non-coding SNP: FxnCode = 13, 15, 53, 55, 6, 73,75; INFO Field = R3, R5, U3, U5, INT, ASS, DSS
coding.classes = dbSNP.df$U3==TRUE|dbSNP.df$U5==TRUE|dbSNP.df$R3==TRUE|dbSNP.df$R5==TRUE|dbSNP.df$INT==TRUE|dbSNP.df$ASS==TRUE|dbSNP.df$DSS==TRUE
snp.filtered = dbSNP.df[coding.classes,]
write.table(snp.filtered[c("CHR", "SNP")], file = sprintf("%s/NWE_data/annotation/noncoding.snp", project.dir), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)

#4. synomymous SNP: FxnCode = 3; INFO Field = SYN
snp.filtered = dbSNP.df[dbSNP.df$SYN == TRUE,]
write.table(snp.filtered[c("CHR", "SNP")], file = sprintf("%s/NWE_data/annotation/synomymous.snp", project.dir), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)

#4. synomymous non-SNP: FxnCode = 41,42,44; INFO Field = NSN, NSM, NSF
snp.filtered = dbSNP.df[dbSNP.df$NSN == TRUE | dbSNP.df$NSM == TRUE| dbSNP.df$NSF == TRUE,]
write.table(snp.filtered[c("CHR", "SNP")], file = sprintf("%s/NWE_data/annotation/nonsynonymous.snp", project.dir), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)

#5. nonsynomymous nonsence: FxnCode = 41; INFO Field = NSN
snp.filtered = dbSNP.df[dbSNP.df$NSN == TRUE,]
write.table(snp.filtered[c("CHR", "SNP")], file = sprintf("%s/NWE_data/annotation/nonsense.snp", project.dir), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)

#6. nonsynomymous missence: FxnCode = 42; INFO Field = NSM
snp.filtered = dbSNP.df[dbSNP.df$NSM == TRUE,]
write.table(snp.filtered[c("CHR", "SNP")], file = sprintf("%s/NWE_data/annotation/missense.snp", project.dir), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)

#7. nonsynomymous non-SNP: FxnCode = 44; INFO Field = NSF
snp.filtered = dbSNP.df[dbSNP.df$NSF == TRUE,] #just 14 SNPs
# write.table(snp.filtered[c("CHR", "SNP")], file = sprintf("%s/NWE_data/annotation/missense.snp", project.dir), sep = "\t", 
#             row.names = FALSE, col.names = TRUE, quote=FALSE)

# ====================================================================
#            Gene annotation for EoH outliers
# ====================================================================
vcf.file = "EoH_outliers_analysis/dbSNP_annotation_filtering.vcf"

#check for info columns  
hdr = scanVcfHeader(vcf.file)
info.header = info(hdr)
info.fields = rownames(info.header)

# info field that contains dbSNP Database codes: R3 R5 SYN NSN NSM NSF U3 U5 INT ASS DSS REF
info.list = c("GENEINFO")
param = ScanVcfParam(info=info.list)

tic("vcf file loading")
vcf = readVcf(vcf.file, param=param)
toc()

vcf.rows = rowRanges(vcf)
df.CHR = data.frame(seqnames(vcf.rows))
colnames(df.CHR) = "CHR"
df.SNP = data.frame(ranges(vcf.rows))
df.SNP = df.SNP[c("names", "start")]
colnames(df.SNP) = c("SNP", "position")

# vcf.fixed = fixed(vcf)
vcf.info = info(vcf)

dbSNP.df = cbind(df.CHR,df.SNP, vcf.info)
head(dbSNP.df)

gene.info = dbSNP.df[!is.na(dbSNP.df$GENEINFO), c("CHR", "GENEINFO")]
gene.info$GENEINFO = sapply(strsplit(gene.info$GENEINFO,":"), `[`, 1)
gene.info = unique(gene.info)
write.table(gene.info, file = sprintf("%s/EoH_outliers_analysis/gene_info.txt", project.dir), sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote=FALSE)
