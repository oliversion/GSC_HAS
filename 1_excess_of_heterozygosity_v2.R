library(config)
config = config::get(file = "config_files/config.yml")

if (config$install.libraries){
  install.packages("ggplot2")
  install.packages("plotly")
  install.packages("nortest")
  install.packages("moments")
  install.packages("tictoc")
  install.packages("GGally")
  install.packages("gridExtra")
}

library(tictoc)

# upload all functions to estimate metrics, plot graphs and test for significance
source(sprintf("%s/0_support_func.R", config$project.dir))

# ====================================================================
#   Create subfolders in the project.dir for saving outputs if needed
# ====================================================================

# structure of project.dir
# - saved_data (to save dataframe with metrics and combined MAF)
# - logs (to save log file with test results and statictics)
# - plots (plots will be sanve to a subfolder corresponding to the "prefix" value)
ifelse(!dir.exists(file.path(config$project.dir, "saved_data")), dir.create(file.path(config$project.dir, "saved_data")), FALSE)
ifelse(!dir.exists(file.path(config$project.dir, "logs")), dir.create(file.path(config$project.dir, "logs")), FALSE)
ifelse(!dir.exists(file.path(config$project.dir, "plots")), dir.create(file.path(config$project.dir, "plots")), FALSE)
plots.path = sprintf("%s/plots", config$project.dir)
ifelse(!dir.exists(file.path(plots.path, config$stat.plots.dir)), dir.create(file.path(plots.path, config$stat.plots.dir)), FALSE)

# ====================================================================
#                               MAF
# ====================================================================

if (!config$recovery.mode){
  # Combined MAF
  # plink --bfile HA_data_1_2_non0MAF --freq --out plink_output/MAF_1_22_non0MAF 
  # MAF for case/controls female/male 
  # plink --bfile HA_data_1_2_non0MAF --freq case-control --filter-males --out plink_output/MAF_cc_male_non0MAF 
  # Combined MAF for female/male
  # plink --bfile HA_data_1_2_non0MAF --freq --filter-males --out plink_output/MAF_male_non0MAF 
  tic("reading plink.MAF data...")
  plink.MAF = read.table(sprintf("%s/%s", config$project.dir, config$plink.MAF.path), header =TRUE, as.is=T)
  toc()
  save(plink.MAF, file = sprintf("%s/saved_data/%sMAF.RData", config$project.dir, config$subgroup.pref))
}

# restore pre-saved data to plink.MAF varialble
tic("reading pre-saved plink.MAF data...")
load(sprintf("%s/saved_data/%sMAF.RData", config$project.dir, config$subgroup.pref))
toc()

#====================================================================
#            Create, save, load the dataframe from plink.hwe
# ====================================================================

dataframe.file = sprintf("%s/saved_data/%sdf.RData", config$project.dir, config$subgroup.pref)
log.file = sprintf("%s/logs/%s%soutput.log", config$project.dir, config$subgroup.pref, config$snp.filter)

if (!config$recovery.mode){
  if (config$is.case.control.analysis){
    
    # plink --bfile NWE_data/HA_data_1_22_non0MAF --hardy  --out plink_output/HWE_cc
    tic("reading plink.hardy data...")
    plink.hardy = read.table(sprintf("%s/%s", config$project.dir, config$plink.hardy.path), header =TRUE, as.is=T)
    toc()
    
    # create dataframe: CHR, SNP, GROUP, GENO, O.HET. , E.HET. , EoH, MAF.bin
    tic("creating dataframe...")
    df = eoh.create.df(data = plink.hardy, filter = plink.MAF)
    toc()
    save(df, file = dataframe.file)
  }
  else {
    #TODO Rawnak
    # add code for creating a dataframe from 7 plink outputs 
    # dataframe structure: CHR, SNP, GROUP, GENO, O.HET. , E.HET. , EoH, MAF.bin
    
    #df = ....
    #save(df, file = dataframe.file)
  }
}

# load pre-saved dataframe: df
# dataframe structure: CHR, SNP, GROUP, GENO, O.HET. , E.HET. , EoH, MAF.bin
tic("reading pre-saved data frame with het metrics...")
load(dataframe.file)
toc()

#====================================================================
#            SNP FILTER
# ====================================================================
if (!config$snp.filter == ''){
   snp.filter = read.table(sprintf("%s", config$snp.filter.path), header =TRUE, as.is=T)
   #filter datadramy by merging
   df = base::merge(x = df, y = snp.filter, by = c("CHR", "SNP"))
 }

#====================================================================
#       EXSESS OF HETEROZYGOSITY and MAF binning tests
# ====================================================================

tic("testing...")
eoh.test(df, log.file, param = config)
toc()

# analysis of the top right ouliers (EoH for HA/controls > 0.12) for EoH metric: MAF, annotation, position
# eoh.outliers.analysis(df)
# eoh.max.dif.analysis(df)
snp.map = read.table(sprintf("%s", config$snp.map.file), header =TRUE, as.is=T)
df.map = base::merge(x = df, y = snp.map, by = c("CHR", "SNP"), all.x = TRUE) #add position to SNP
df.map.HA = df.map[df.map$GROUP == "HA",]
df.map.controls = df.map[df.map$GROUP == "controls",]

png("manhattanEoH_HA.png", width=950, height=500)
print(manhattan.plot(df.map.HA$CHR, df.map.HA$position, df.map.HA$EoH, ylab = "EoH: HA"))
dev.off()

png("manhattanEoH_controls.png", width=950, height=500)
print(manhattan.plot(df.map.controls$CHR, df.map.controls$position, 
                     df.map.controls$EoH, ylab = "EoH: controls"))
dev.off()

df.map.HA$dif = df.map.HA$EoH-df.map.controls$EoH
png("manhattanEoH_dif.png", width=950, height=500)
print(manhattan.plot(df.map.HA$CHR, df.map.HA$position, 
                     df.map.HA$dif, ylab = "EoH: HA-controls", 
                     col=c("blue4", "orange3")))
dev.off()

# Chromosome 6
df.map.HA.6 = df.map.HA[df.map.HA$CHR==6,]
ann = rep(1, length(df.map.HA.6$dif))
ann[df.map.HA.6$CHR == 6 & df.map.HA.6$position>=28477797 & df.map.HA.6$position<=33448354]= 2
ann = factor(ann, levels=1:2, labels=c("","HLA"))

png("manhattanEoH_dif_6.png", width=950, height=500)
print(manhattan.plot(df.map.HA.6$CHR, df.map.HA.6$position, 
                     df.map.HA.6$dif, ylab = "EoH for 6 CHR: HA-controls", 
                     annotate=ann, sig.level=0.4, col=c("blue4", "orange3")))
dev.off()

#HLA
df.map.HLA = df.map.HA[df.map.HA$CHR == 6 & df.map.HA$position>=28477797 & df.map.HA$position<=33448354,]
png("manhattanEoH_dif_HLA.png", width=950, height=500)
print(manhattan.plot(df.map.HLA$CHR, df.map.HLA$position, 
                     df.map.HLA$dif, ylab = "EoH for HLA: HA-controls", sig.level=0.4,
                     col=c("blue4", "orange3")))
dev.off()

#draw plot with annotation
ann = rep(1, length(df.map.HLA$dif))
ann[with(df.map.HLA, EoH >=0.4 & EoH<= -0.4)]= 1
ann = factor(ann, levels=1:2, labels=c("","outliers"))
png("manhattanEoH_dif_HLA_outliers.png", width=950, height=500)
print(manhattan.plot(df.map.HLA$CHR, df.map.HLA$position, 
                     df.map.HLA$dif, annotate=ann, ylab = "EoH for HLA: HA-controls", sig.level=0.4,
                     col=c("blue4", "orange3")))
dev.off()

#====================================================================
#       EXSESS OF HOMOZYGOSITY FOR MINOR AND MAJOR ALLELE
# ====================================================================

if (config$do.homozygosity.test){
  # http://www.uwyo.edu/dbmcd/popecol/maylects/fst.html
  log.file = sprintf("%s/logs/%s%sHOMOZ_output.log", config$project.dir, config$subgroup.pref, config$snp.filter)
  
  plots.path = sprintf("%s/plots/%s", config$project.dir, config$stat.plots.dir)
  ifelse(!dir.exists(file.path(plots.path, "HOMOZ")), dir.create(file.path(plots.path, "HOMOZ")), FALSE)
  
  tic("testing...")
  df = homoz.test(df, log.file, param = config)
  toc()
}

#====================================================================
#      LINE PLOT and HEAT MAP and smoothing per chromosome
# ====================================================================
if (config$save.plots.chr){
  plot.chr(df, param = config, varname = "EoH")
  plot.chr(df, param = config, varname = "EoHOMOZ")
  plot.chr(df, param = config, varname = "EoHOMOZ.major")
}

#====================================================================
#      ODD RATION ANALYSIS
# ====================================================================
if (config$do.odd.ratio.test){
  df.odds.ratio = odd.ratio.test(df, param = config)
  
}

# in case of any error all counnections (opened with sink) should be closed otherwise output will be written in an opened log file 
# closeAllConnections() 