library(plyr)
library(ggplot2)
library(gtable)

# ====================================================================
# 1.          Excess of heterozygosity 
# ====================================================================
eoh.create.df = function(data, filter){
  #check for zero MAF
  printf("#zeros for Het expected cases: %i\n", sum(data$E.HET.[data $TEST=='AFF']==0))
  printf("#zeros for Het expected controls: %i\n", sum(data$E.HET.[data $TEST=='UNAFF']==0))
  
  df = data[data$TEST != 'ALL', c("CHR", "SNP", "TEST", "GENO", "O.HET." , "E.HET.")]
  colnames(df) = c("CHR", "SNP", "GROUP", "GENO", "O.HET." , "E.HET." )
  df$GROUP[df$GROUP== "AFF"] = "HA"
  df$GROUP[df$GROUP== "UNAFF"] = "controls"
  df$EoH =  (df$O.HET. - df$E.HET.)/df$E.HET.
  
  df$HOM_minor = as.numeric(sapply(strsplit(df$GENO,"/"), `[`, 1))
  df$HET       = as.numeric(sapply(strsplit(df$GENO,"/"), `[`, 2))
  df$HOM_major = as.numeric(sapply(strsplit(df$GENO,"/"), `[`, 3))
  df$total     = df$HOM_minor + df$HET + df$HOM_major
  
  #for a minor allele
  df$E.HOMOZ   = ((2*df$HOM_minor + df$HET)/(2*df$total))^2
  df$O.HOMOZ   = df$HOM_minor/df$total
  df$EoHOMOZ   = (df$O.HOMOZ - df$E.HOMOZ)/df$E.HOMOZ
  
  #for a major allele
  df$E.HOMOZ.major = 1-df$E.HOMOZ-df$E.HET.
  df$O.HOMOZ.major = df$HOM_major/df$total
  df$EoHOMOZ.major = (df$O.HOMOZ.major - df$E.HOMOZ.major)/df$E.HOMOZ.major
  
  # add the info about MAF bins
  maf.bins.df = get.maf.bins(plink.MAF)
  
  # merging two dataframes on CHR and SNP
  df = base::merge(x = df, y = maf.bins.df, by = c("CHR", "SNP"), all.x = TRUE)
  return(df)
}
eoh.test = function(in.df, log.file, param){
  # colnames(df): "CHR" "SNP" "GROUP" "GENO" "O.HET." "E.HET." "EoH" "MAF.bin"
  eoh.summary = data_summary(in.df, varname="EoH", groupnames=c("GROUP", "MAF.bin"))
  eoh.summary$MAF.bin=as.factor(eoh.summary$MAF.bin)
  
  het.obs.summary = data_summary(in.df, varname="O.HET.", groupnames=c("GROUP", "MAF.bin"))
  het.obs.summary$MAF.bin=as.factor(het.obs.summary$MAF.bin)
  
  # Start writing to the log file
  sink(log.file)
  
  cat("=======================================================\n")
  cat("EXPERIMENTAL SET-UP\n")
  cat("=======================================================\n")
  cat(paste(names(param), param, sep = ":", collapse = "\n"))
  printf("\nTOTAL # of SNPs: %i\n", length(unique(in.df$SNP)))
  
  cat("\n=======================================================\n")
  cat("OBSERVED HETEROZYGOSITY\n")
  cat("=======================================================\n")

  show.df.summary(het.obs.summary)
    
  cat("=======================================================\n")
  cat("EXCESS OF HETEROZYGOSITY\n")
  cat("=======================================================\n")

  show.df.summary(eoh.summary)
  
  if (param$do.stat.test){
    # The unpaired two-samples Wilcoxon test (also known as Wilcoxon rank sum test or Mann-Whitney test) 
    # is a non-parametric alternative to the unpaired two-samples t-test, which can be used to compare two independent groups 
    # of samples. It’s used when the data are not normally distributed.
    if (param$is.case.control.analysis){
      print(wilcox.test(in.df$EoH[in.df$GROUP == "controls"], in.df$EoH[in.df$GROUP == "HA"], exact = FALSE, alternative = "less"))
      print(fligner.test(in.df$EoH[in.df$GROUP == "controls"], in.df$EoH[in.df$GROUP == "HA"], alternative = "less"))
    }
    else{
      #TODO RAWNAK check for multiple groups
      print(wilcox.test(EoH~GROUP, in.df, exact = FALSE))
      print(fligner.test(EoH~GROUP, in.df, exact = FALSE)) 
    }
  }
  
  if (param$save.plots == TRUE) {
    if (param$is.case.control.analysis){
      plot.df.stats(in.df, param)
    }
    
    labels = get.maf.bins.names()
    maf.bin.count = table(in.df$MAF.bin)/2
    for (i in 1:length(labels)){
      labels[i] = sprintf("%s\n%s SNPs", labels[i], maf.bin.count[i])
    }
    plot.maf.binning.stat(het.obs.summary, varname="HET_Obs", param, labels)
    plot.maf.binning.stat(eoh.summary, varname="EoH", param, labels)
    plot.density(in.df, eoh.summary, param, varname="EoH")
    plot.eoh(in.df, param, varname="EoH")
  }
  
  cat("=======================================================\n")
  cat("MAF binning\n")
  cat("=======================================================\n")

  maf.bins.name = get.maf.bins.names()
  for (i in 1:length(maf.bins.name)){
    maf.bin.df = in.df[in.df$MAF.bin == i,]
    
    printf("\nMAF bin: %s. # of SNPs: %i\n", maf.bins.name[i], length(unique(maf.bin.df$SNP)))
    cat("=======================================================\n")
    
    if (param$do.stat.test) {
      if (param$is.case.control.analysis){
        print(wilcox.test(maf.bin.df$EoH[maf.bin.df$GROUP == "controls"], maf.bin.df$EoH[maf.bin.df$GROUP == "HA"], exact = FALSE, alternative = "less"))
        print(fligner.test(maf.bin.df$EoH[maf.bin.df$GROUP == "controls"], maf.bin.df$EoH[maf.bin.df$GROUP == "HA"], alternative = "less"))
      }
      else{
        #TODO RAWNAK check for multiple groups
        print(wilcox.test(EoH~GROUP, maf.bin.df, exact = FALSE))
        print(fligner.test(EoH~GROUP, maf.bin.df, exact = FALSE)) 
      }
    }
    
    if (param$save.plots == TRUE) {
      maf.bin.df.summary = data_summary(maf.bin.df, varname="EoH", groupnames=c("GROUP"))
      plot.density(maf.bin.df, maf.bin.df.summary, param, varname="EoH", maf.bin = i)
      plot.eoh(maf.bin.df, param, maf.bin = i)
    }
  }
  # Stop writing to the log file
  sink()
}
homoz.test = function(in.df, log.file, param){
  # colnames(in.df): "CHR" "SNP" "GROUP" "GENO" "O.HET." "E.HET." "EoH" "MAF.bin"

  #for a minor allele
  EoHOMOZ.summary = data_summary(in.df, varname="EoHOMOZ", groupnames=c("GROUP", "MAF.bin"))
  EoHOMOZ.summary$MAF.bin=as.factor(EoHOMOZ.summary$MAF.bin)
  
  O.HOMOZ.summary = data_summary(in.df, varname="O.HOMOZ", groupnames=c("GROUP", "MAF.bin"))
  O.HOMOZ.summary$MAF.bin=as.factor(O.HOMOZ.summary$MAF.bin)
  
  #for a major allele
  EoHOMOZ.major.summary = data_summary(in.df, varname="EoHOMOZ.major", groupnames=c("GROUP", "MAF.bin"))
  EoHOMOZ.major.summary$MAF.bin=as.factor(EoHOMOZ.major.summary$MAF.bin)
  
  # Start writing to the log file
  sink(log.file)
  cat("=======================================================\n")
  cat("EXPERIMENTAL SET-UP\n")
  cat("=======================================================\n")
  cat(paste(names(param), param, sep = ":", collapse = "\n"))
  printf("\nTOTAL # of SNPs: %i\n", length(unique(in.df$SNP)))
  
  cat("\n=======================================================\n")
  cat("OBSERVED HOMOZYGOSITY FOR MINOR ALLELE\n")
  cat("=======================================================\n")
  
  show.df.summary(O.HOMOZ.summary)
  
  cat("\n=======================================================\n")
  cat("EXCESS OF HOMOZYGOSITY FOR MINOR ALLELE\n")
  cat("=======================================================\n")
  
  show.df.summary(EoHOMOZ.summary)
  
  cat("\n=======================================================\n")
  cat("EXCESS OF HOMOZYGOSITY FOR MAJOR ALLELE\n")
  cat("=======================================================\n")
  
  show.df.summary(EoHOMOZ.major.summary)
  
  if (param$do.stat.test){
    # The unpaired two-samples Wilcoxon test (also known as Wilcoxon rank sum test or Mann-Whitney test) 
    # is a non-parametric alternative to the unpaired two-samples t-test, which can be used to compare two independent groups 
    # of samples. It’s used when the data are not normally distributed.
    if (param$is.case.control.analysis){
      print(wilcox.test(in.df$EoHOMOZ[in.df$GROUP == "controls"], in.df$EoHOMOZ[in.df$GROUP == "HA"], exact = FALSE, alternative = "less"))
      print(fligner.test(in.df$EoHOMOZ[in.df$GROUP == "controls"], in.df$EoHOMOZ[in.df$GROUP == "HA"], alternative = "less"))
    }
    else{
      #TODO RAWNAK check for multiple groups
      print(wilcox.test(EoHOMOZ~GROUP, in.df, exact = FALSE))
      print(fligner.test(EoHOMOZ~GROUP, in.df, exact = FALSE)) 
    }
  }
  
  if (param$save.plots == TRUE) {
    labels = get.maf.bins.names()
    maf.bin.count = table(in.df$MAF.bin)/2
    for (i in 1:length(labels)){
      labels[i] = sprintf("%s\n%s SNPs", labels[i], maf.bin.count[i])
    }
    plot.maf.binning.stat(O.HOMOZ.summary, varname="HOMOZ_Obs", param, labels)
    plot.maf.binning.stat(EoHOMOZ.summary, varname="EoHOMOZ", param, labels)
    plot.maf.binning.stat(EoHOMOZ.major.summary, varname="EoHOMOZ.major", param, labels)
    plot.density(in.df, EoHOMOZ.summary, param, varname="EoHOMOZ")
    plot.density(in.df, EoHOMOZ.major.summary, param, varname="EoHOMOZ.major")
    plot.eoh(in.df, param, varname="EoHOMOZ")
    plot.eoh(in.df, param, varname="EoHOMOZ.major")
  }
  
  cat("=======================================================\n")
  cat("MAF binning\n")
  cat("=======================================================\n")
  
  maf.bins.name = get.maf.bins.names()
  for (i in 1:3){
    maf.bin.df = in.df[in.df$MAF.bin == i,]
    
    printf("\nMAF bin: %s. # of SNPs: %i\n", maf.bins.name[i], length(unique(maf.bin.df$SNP)))
    cat("=======================================================\n")
    
    if (param$do.stat.test) {
      if (param$is.case.control.analysis){
        print(wilcox.test(maf.bin.df$EoHOMOZ[maf.bin.df$GROUP == "controls"], maf.bin.df$EoHOMOZ[maf.bin.df$GROUP == "HA"], exact = FALSE, alternative = "less"))
        print(fligner.test(maf.bin.df$EoHOMOZ[maf.bin.df$GROUP == "controls"], maf.bin.df$EoHOMOZ[maf.bin.df$GROUP == "HA"], alternative = "less"))
      }
      else{
        #TODO RAWNAK check for multiple groups
        print(wilcox.test(EoHOMOZ~GROUP, maf.bin.df, exact = FALSE))
        print(fligner.test(EoHOMOZ~GROUP, maf.bin.df, exact = FALSE)) 
      }
    }
    
    if (param$save.plots == TRUE) {
      plot.density(maf.bin.df, data_summary(maf.bin.df, varname="EoHOMOZ", groupnames=c("GROUP")), param, 
                   varname="EoHOMOZ", maf.bin = i)
      plot.density(maf.bin.df, data_summary(maf.bin.df, varname="EoHOMOZ.major", groupnames=c("GROUP")), param, 
                   varname="EoHOMOZ.major", maf.bin = i)
      plot.eoh(maf.bin.df, param, varname="EoHOMOZ", maf.bin = i)
      plot.eoh(maf.bin.df, param, varname="EoHOMOZ.major", maf.bin = i)
    }
  }
  # Stop writing to the log file
  sink()
  return(in.df)
}
odd.ratio.test = function(in.df, param){
  in.df$R.AA     = in.df$HOM_major/in.df$HOM_minor
  in.df$R.Aa     = in.df$HET/in.df$HOM_minor
  in.df = base::merge(x = in.df, y = plink.MAF, by = c("CHR", "SNP"), all.x = TRUE)
  
  df.odds.ratio = in.df[in.df$GROUP=="HA",]
  df.odds.ratio$R.AA.controls = in.df$R.AA[in.df$GROUP=="controls"]
  df.odds.ratio$R.Aa.controls = in.df$R.Aa[in.df$GROUP=="controls"]
  df.odds.ratio$OR.AA = df.odds.ratio$R.AA/df.odds.ratio$R.AA.controls
  df.odds.ratio$OR.Aa = df.odds.ratio$R.Aa/df.odds.ratio$R.Aa.controls
  return(df.odds.ratio)
}

# TODO Rawnak this function should be adaped if needed. Used onle for a report
plot.df.stats = function(in.df, param){
  # plot statistics for the df: heterozygosity observes and expected, call rate, histograms etc. 
  # colnames(df): "CHR" "SNP" "GROUP" "GENO" "O.HET." "E.HET." "EoH" "MAF.bin"
  plots.path = sprintf("%s/plots", param$project.dir)
  file.path = sprintf("%s/%s/", plots.path, param$stat.plots.dir)
  
  df.HA       = in.df[in.df$GROUP == "HA", ]
  df.controls = in.df[in.df$GROUP == "controls", ]

  png(sprintf("%s/plots/%s/%s%sHET_obs_cc.png", param$project.dir, param$stat.plots.dir, param$subgroup.pref, param$snp.filter), 
              width = 350, height = 350)
  plot(df.HA$O.HET., df.controls$O.HET., main = "HET_obs HA vs controls", xlab = "HA", ylab = "controls", pch='.')
  dev.off() 
  
  png(sprintf("%s/plots/%s/%s%sHET_exp_cc.png", param$project.dir, param$stat.plots.dir, param$subgroup.pref, param$snp.filter),
              width = 350, height = 350)
  plot(df.HA$E.HET., df.controls$E.HET., main = "HET_exp HA vs controls", xlab = "HA", ylab = "controls", pch='.')
  dev.off()
  
  png(sprintf("%s/plots/%s/%s%sEoH_cc.png", param$project.dir, param$stat.plots.dir, param$subgroup.pref, param$snp.filter),
              width = 350, height = 350)
  plot(df.HA$EoH, df.controls$EoH, main = "EoH HA vs controls", xlab = "HA", ylab = "controls", pch='.')
  dev.off()
  
  png(sprintf("%s/plots/%s/%s%sEoHOMOZ_cc.png", param$project.dir, param$stat.plots.dir, param$subgroup.pref, param$snp.filter),
      width = 350, height = 350)
  plot(df.HA$EoHOMOZ, df.controls$EoHOMOZ, main = "EoHOMOZ HA vs controls", xlab = "HA", ylab = "controls", pch='.')
  dev.off()
  
  png(sprintf("%s/plots/%s/%s%sEoHOMOZ.major_cc.png", param$project.dir, param$stat.plots.dir, param$subgroup.pref, param$snp.filter),
      width = 350, height = 350)
  plot(df.HA$EoHOMOZ.major, df.controls$EoHOMOZ.major, main = "EoHOMOZ.major HA vs controls", xlab = "HA", ylab = "controls", pch='.')
  dev.off()
}
plot.maf.binning.stat = function(in.df.summary, varname, param, labels=NA) {
  file.path = sprintf("%s/plots/%s/", param$project.dir, param$stat.plots.dir)
  if (varname == "HOMOZ_Obs" | varname =="EoHOMOZ" | varname =="EoHOMOZ.major"){
    file.path = sprintf("%s/%s/", plots.path, "HOMOZ")
  }
  
  # Line plot for MEAN and standart error
  p = ggplot(in.df.summary, aes(x=MAF.bin, y=mean, group=GROUP, color=GROUP)) +
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(0.05)) +
    scale_fill_brewer(palette="Paired") + theme_minimal() +
    theme(axis.text=element_text(size=5), axis.title=element_text(size=6)) + 
    labs(title = sprintf("Mean of %s per MAF bin", varname), x="MAF bin", y = varname) +
    scale_x_discrete(labels = labels)
  ggsave(sprintf("%s%s%s_mafbinned_mean_se.pdf", file.path, param$subgroup.pref, varname))

  # Line plot for MEDIAN and standart error
  p = ggplot(in.df.summary, aes(x=MAF.bin, y=median, group=GROUP, color=GROUP)) +
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=median-se, ymax=median+se), width=.2, position=position_dodge(0.05)) +
    scale_fill_brewer(palette="Paired") + theme_minimal() +
    theme(axis.text=element_text(size=5), axis.title=element_text(size=6)) + 
    labs(title = sprintf("Median of %s per MAF bin", varname), x="MAF bin", y = varname)+
    scale_x_discrete(labels = labels)
  ggsave(sprintf("%s%s%s_mafbinned_median_se.pdf", file.path, param$subgroup.pref, varname))
  
  # # Line plot for MEAN and standart deviation
  p = ggplot(in.df.summary, aes(x=MAF.bin, y=mean, group=GROUP, color=GROUP)) +
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.05)) +
    scale_fill_brewer(palette="Paired") + theme_minimal() +
    theme(axis.text=element_text(size=5), axis.title=element_text(size=6)) + 
    labs(title = sprintf("Mean of %s per MAF bin", varname), x="MAF bin", y = varname) +
    scale_x_discrete(labels = labels)
  ggsave(sprintf("%s%s%s_mafbinned_mean_sd.pdf", file.path, param$subgroup.pref, varname))

  # # Line plot for MEDIAN and standart deviation
  p = ggplot(in.df.summary, aes(x=MAF.bin, y=median, group=GROUP, color=GROUP)) +
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.2, position=position_dodge(0.05)) +
    scale_fill_brewer(palette="Paired") + theme_minimal() +
    theme(axis.text=element_text(size=5), axis.title=element_text(size=6)) + 
    labs(title = sprintf("Median of %s per MAF bin", varname), x="MAF bin", y = varname)+
    scale_x_discrete(labels = labels)
  ggsave(sprintf("%s%s%s_mafbinned_median_sd.pdf", file.path, param$subgroup.pref, varname))
 
  # barplot for MEAN and standart error
  p <- ggplot(in.df.summary, aes(x=MAF.bin, y=mean, fill=GROUP, label = mean)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=sprintf("%.4f",mean)), vjust=-1, size = 2, position = position_dodge(width = 1), hjust = -0.1) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + 
    scale_fill_brewer(palette="Paired") + theme_minimal() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=6)) + 
    labs(title = sprintf("Mean of %s per MAF bin", varname), x="MAF bin", y = varname) +
    scale_x_discrete(labels = labels)
  ggsave(sprintf("%s%s%s_mafbinned_mean_se_barplot.pdf", file.path, param$subgroup.pref, varname))
  
  # barplot for MEDIAN and standart error
  p <- ggplot(in.df.summary, aes(x=MAF.bin, y=median, fill=GROUP, label = median)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=sprintf("%.4f",median)), vjust=-1, size = 2, position = position_dodge(width = 1), hjust = -0.1) +
    geom_errorbar(aes(ymin=median-se, ymax=median+se), width=.2, position=position_dodge(.9)) +
    scale_fill_brewer(palette="Paired") + theme_minimal() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=6)) + 
    labs(title = sprintf("Median of %s per MAF bin", varname), x="MAF bin", y = varname) +
    scale_x_discrete(labels = labels)
  ggsave(sprintf("%s%s%s_mafbinned_median_se_barplot.pdf", file.path, param$subgroup.pref, varname))
}
plot.eoh = function(in.df, param, varname="EoH", maf.bin = 0){
  
  plots.path = sprintf("%s/plots", param$project.dir)
  file.path = sprintf("%s/%s/", plots.path, param$stat.plots.dir)
  if (varname == "HOMOZ_Obs" | varname =="EoHOMOZ" |  varname =="EoHOMOZ.major"){
    file.path = sprintf("%s/%s/", file.path, "HOMOZ")
  }
  
  file.prefix = sprintf("%s%s_%s", param$subgroup.pref, varname, param$snp.filter)
  if (maf.bin > 0){
    file.prefix = paste0(file.prefix,"_mafbin", maf.bin, collapse = NULL)
  }
  
  # Boxplots for each group
  p = ggplot(in.df, aes_string(x="GROUP", y=varname, color="GROUP")) + geom_boxplot()
  ggsave(sprintf("%s%sboxplot.png", file.path, file.prefix))
  
  }

plot.density = function(in.df, in.df.summary, param, varname="EoH", maf.bin = 0){
  file.path = sprintf("%s/plots/%s/", param$project.dir, param$stat.plots.dir)
  if (varname == "HOMOZ_Obs" | varname =="EoHOMOZ" |  varname =="EoHOMOZ.major"){
    file.path = sprintf("%s/%s/", plots.path, "HOMOZ")
  }
  
  file.prefix = sprintf("%s%s_%s", param$subgroup.pref, varname, param$snp.filter)
  if (maf.bin > 0){
    file.prefix = paste0(file.prefix,"_mafbin", maf.bin, collapse = NULL)
  }
  #1
  ggplot(in.df, aes_string(x=varname, color="GROUP")) + xlab(varname) +
    geom_density() +
    geom_vline(data=in.df.summary, aes(xintercept=mean,  colour=GROUP), linetype="dashed", size=0.5)
  ggsave(sprintf("%s%s_geom_dens.pdf", file.path, file.prefix))

  #2
  ggplot(in.df, aes_string(x=varname, fill="GROUP")) + geom_density(alpha=.3) + xlab(varname) +
    geom_vline(data=in.df.summary, aes(xintercept=mean,  colour=GROUP), linetype="dashed", size=0.5)
  ggsave(sprintf("%s%s_geom_dens_color.pdf", file.path, file.prefix))

  #3
  ggplot(in.df, aes_string(x=varname, color="GROUP")) + xlab(varname) +
    geom_histogram(binwidth=.001, alpha=.5, position="identity") +
    geom_vline(data=in.df.summary, aes(xintercept=mean,  colour=GROUP), linetype="dashed", size=0.5)
  ggsave(sprintf("%s%s_geom_hist.pdf", file.path, file.prefix))
}
plot.chr = function(in.df, param, varname = "EoH"){
  snp.map = read.table(sprintf("%s", param$snp.map.file), header =TRUE, as.is=T)
  df.map = base::merge(x = in.df, y = snp.map, by = c("CHR", "SNP"), all.x = TRUE) #add position to SNP
  chr.list = unique(df.map$CHR)
  
  file.path = sprintf("%s/plots/%s/", param$project.dir, param$stat.plots.dir)
  if (varname == "HOMOZ_Obs" | varname =="EoHOMOZ" |  varname =="EoHOMOZ.major"){
    file.path = sprintf("%s/%s/", plots.path, "HOMOZ")
  }
  
  # plot EoH separately each chromosome for HA, controls and difference
  for (i in chr.list){
    file.prefix = sprintf("%s%s_%s", param$subgroup.pref, varname, param$snp.filter)
    file.prefix = paste0(file.prefix,"_chr", i, collapse = NULL)
   
    chr.data = df.map[df.map$GROUP=="HA" & df.map$CHR==i, c("CHR", varname, "position")]
    colnames(chr.data) = c("CHR", "EoH", "position")
    chr.data$EoH.controls =  df.map[df.map$GROUP=="controls" & df.map$CHR==i, varname]
    chr.data$EoH.dif = chr.data$EoH - chr.data$EoH.controls 
    chr.data = chr.data[with(chr.data, order(position)),]
    
    #http://r-statistics.co/Loess-Regression-With-R.html
    chr.data$smoothed.EoH.HA  = predict(loess(EoH ~ position, data=chr.data, span=param$span)) 
    chr.data$smoothed.EoH.controls = predict(loess(EoH.controls ~ position, data=chr.data, span=param$span)) 
    chr.data$smoothed.EoH.dif = predict(loess(EoH.dif ~ position, data=chr.data, span=param$span)) 
    
    options("scipen"=100, "digits"=4)
    g1 = ggplot(chr.data, aes(position)) +
      geom_line(aes(y = EoH, colour = "EoH for HA"), size = 0.1) +
      geom_line(aes(y = smoothed.EoH.HA, colour = "EoH trend for HA"), size = 0.1) +
      scale_colour_manual(values=c("gray86", "red")) + 
      guides(shape = guide_legend(override.aes = list(size = 0.4)),color = guide_legend(override.aes = list(size = 0.4))) +
      theme(axis.text=element_text(size=5), axis.title=element_text(size=6)) + 
      geom_hline(yintercept = 0, color="black", size = 0.1)
      #+geom_vline(xintercept=25500000, color="black", linetype="dashed", size=0.1) +
      #geom_vline(xintercept=33500000, color="black", linetype="dashed", size=0.1)
    
    g2 = ggplot(chr.data, aes(position)) +
      geom_line(aes(y = EoH.controls, colour = "EoH for controls"), size = 0.1) +
      geom_line(aes(y = smoothed.EoH.controls, colour = "EoH trend for controls"), size = 0.1) +
      scale_colour_manual(values=c("gray86", "red")) + 
      guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5))) +
      theme(axis.text=element_text(size=5), axis.title=element_text(size=6)) + 
      geom_hline(yintercept = 0, color="black", size = 0.1)
      
    g3 = ggplot(chr.data, aes(position)) +
      geom_line(aes(y = EoH.dif, colour = "EoH:HA-controls"), size = 0.1) +
      geom_line(aes(y = smoothed.EoH.dif, colour = "EoH:HA-controls trend"), size = 0.1) +
      scale_colour_manual(values=c("gray86", "red")) + 
      guides(shape = guide_legend(override.aes = list(size = 0.5)),color = guide_legend(override.aes = list(size = 0.5))) +
      theme(axis.text=element_text(size=5), axis.title=element_text(size=6)) + 
      geom_hline(yintercept = 0, color="black", size = 0.1)
    ggsave(sprintf("%s%s_line.png", file.path, file.prefix), arrangeGrob(g1, g2, g3))
    
    ggplot(chr.data, aes(position)) +
      theme(axis.text=element_text(size=5), axis.title=element_text(size=6)) + 
      geom_line(aes(y = smoothed.EoH.HA, colour = "EoH trend for HA"), size = 0.1) +
      geom_line(aes(y = smoothed.EoH.controls, colour = "EoH trend for controls"), size = 0.1) +
      geom_line(aes(y = smoothed.EoH.dif, colour ="EoH:HA-controls trend"), size = 0.3)
    ggsave(sprintf("%s%s_line_all.png", file.path, file.prefix))
    
    #https://www.biostars.org/p/336999/
    g1 = ggplot(data=chr.data, aes(x=position, y=1)) +
      facet_grid(CHR ~ ., switch='y') +
      geom_tile(aes(fill=smoothed.EoH.HA))  +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.text.y = element_text(angle = 180)) +
      guides(fill=guide_legend(title=sprintf("%s for HA", varname)), shape = guide_legend(override.aes = list(size = 0.4)),
             color = guide_legend(override.aes = list(size = 0.5))) +
      theme(legend.title = element_text(size = 5), legend.text = element_text(size = 5)) +
      theme(axis.text=element_text(size=5), axis.title=element_text(size=6)) + 
      scale_fill_gradientn(colours= rainbow(3),
                           breaks=seq(min(chr.data$smoothed.EoH.HA),max(chr.data$smoothed.EoH.HA),
                                      (max(chr.data$smoothed.EoH.HA)-min(chr.data$smoothed.EoH.HA))/4))
    
    g2 = ggplot(data=chr.data, aes(x=position, y=1)) +
      facet_grid(CHR ~ ., switch='y') +
      geom_tile(aes(fill=smoothed.EoH.controls))  +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.text.y = element_text(angle = 180)) +
      guides(fill=guide_legend(title=sprintf("%s for controls", varname)), shape = guide_legend(override.aes = list(size = 0.4)),
             color = guide_legend(override.aes = list(size = 0.5))) +
      theme(legend.title = element_text(size = 5), legend.text = element_text(size = 5)) +
      theme(axis.text=element_text(size=5), axis.title=element_text(size=6)) + 
      scale_fill_gradientn(colours= rainbow(3),
                           breaks=seq(min(chr.data$smoothed.EoH.controls),max(chr.data$smoothed.EoH.controls),
                                      (max(chr.data$smoothed.EoH.controls)-min(chr.data$smoothed.EoH.controls))/4))
    
    g3 = ggplot(data=chr.data, aes(x=position, y=1)) +
      facet_grid(CHR ~ ., switch='y') +
      geom_tile(aes(fill=smoothed.EoH.dif))  +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.text.y = element_text(angle = 180)) +
      guides(fill=guide_legend(title=sprintf("%s:HA-controls", varname)), shape = guide_legend(override.aes = list(size = 0.4)),
             color = guide_legend(override.aes = list(size = 0.5))) +
      theme(legend.title = element_text(size = 5), legend.text = element_text(size = 5)) +
      theme(axis.text=element_text(size=5), axis.title=element_text(size=6)) + 
      scale_fill_gradientn(colours= rainbow(3),
                           breaks=seq(min(chr.data$smoothed.EoH.dif),max(chr.data$smoothed.EoH.dif),
                                      (max(chr.data$smoothed.EoH.dif)-min(chr.data$smoothed.EoH.dif))/4))
    ggsave(sprintf("%s%s_heatmap.pdf", file.path, file.prefix), arrangeGrob(g1, g2, g3))
  }
}

# ====================================================================
#            UTILITIES
# ====================================================================
printf = function(...) cat(sprintf(...))
get.maf.bins = function(in.filter){
  in.filter$MAF.bin = 0
  in.filter$MAF.bin[in.filter$MAF< 0.01] = 1
  in.filter$MAF.bin[(in.filter$MAF< 0.05)&(in.filter$MAF >= 0.01)] = 2
  in.filter$MAF.bin[(in.filter$MAF< 0.1) &(in.filter$MAF >= 0.05)] = 3
  in.filter$MAF.bin[(in.filter$MAF< 0.25)&(in.filter$MAF >= 0.1)] = 4
  in.filter$MAF.bin[(in.filter$MAF<= 0.5) &(in.filter$MAF >= 0.25)] = 5
  
  in.filter$MAF.bin=as.factor(in.filter$MAF.bin)
  return(in.filter[c("CHR", "SNP", "MAF.bin")])
}
get.maf.bins.names = function(){
  maf.bins.name = list("[0%,1%)","[1%,5%)","[5%,10%)", "[10%,25%)", "[25%,50%]")
}
get.non0.maf = function(plink.maf.cc, plink.MAF.cc.female, plink.MAF.cc.male) {
  
  maf.0.cases = plink.maf.cc$MAF_A == 0
  printf("MAF== 0 for HA: %i\n", sum(maf.0.cases))
  
  maf.0.controls = plink.maf.cc$MAF_U == 0
  printf("MAF== 0 for controls: %i\n", sum(maf.0.controls))
  
  maf.0.cases.f = plink.MAF.cc.female$MAF_A == 0
  printf("MAF== 0 for female HA: %i\n", sum(maf.0.cases.f))
  
  maf.0.controls.f = plink.MAF.cc.female$MAF_U == 0
  printf("MAF== 0 for female controls: %i\n", sum(maf.0.controls.f))
  
  maf.0.cases.m = plink.MAF.cc.male$MAF_A == 0
  printf("MAF== 0 for male HA: %i\n", sum(maf.0.cases.m))
  
  maf.0.controls.m = plink.MAF.cc.male$MAF_U == 0
  printf("MAF== 0 for male controls: %i\n", sum(maf.0.controls.m))
  
  snp.use = (plink.MAF.cc.male$MAF_U == 0)|(plink.MAF.cc.male$MAF_A == 0)|(plink.MAF.cc.female$MAF_A == 0)|(plink.MAF.cc.female$MAF_U == 0)
  snp.use = !snp.use
}

# Calculate the mean and the standard deviation for each group
data_summary = function(data, varname, groupnames){
  summary_func = function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE), sd = sd(x[[col]], na.rm=TRUE), 
             median = median(x[[col]], na.rm=TRUE), se = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
  }
  data_sum = ddply(data, groupnames, .fun=summary_func, varname)
  # data_sum = rename(data_sum, c("mean" = varname))
  return(data_sum)
}
show.df.summary = function(df.summary,  varname){
  df.to.show = reshape(df.summary[c("GROUP", "MAF.bin", "mean")], direction = "wide", idvar = "GROUP", timevar = "MAF.bin")
  rownames(df.to.show) <- df.to.show$GROUP
  df.to.show = df.to.show[,2:dim(df.to.show)[2]]
  
  #TODO RAWNAK check how difference works with multiple groups
  diff.df = data.frame(diff(as.matrix(df.to.show)))
  
  colnames(df.to.show) = unlist(get.maf.bins.names())
  colnames(diff.df)    = unlist(get.maf.bins.names())
  
  printf("\n MEAN: \n")
  cat("=======================================================\n")
  print(df.to.show)
  
  printf("\n MEAN DIFFERENCE: (HA - controls): \n")
  print(diff.df)
  
  # median summary
  df.to.show = reshape(df.summary[c("GROUP", "MAF.bin", "median")], direction = "wide", idvar = "GROUP", timevar = "MAF.bin")
  rownames(df.to.show) <- df.to.show$GROUP
  df.to.show = df.to.show[,2:dim(df.to.show)[2]]
  
  #TODO RAWNAK check how difference works with multiple groups
  diff.df = data.frame(diff(as.matrix(df.to.show)))
  
  colnames(df.to.show) = unlist(get.maf.bins.names())
  colnames(diff.df)    = unlist(get.maf.bins.names())
  printf("\n MEDIAN: \n")
  cat("=======================================================\n")
  print(df.to.show)
  
  printf("\n MEDIAN DIFFERENCE: (HA - controls): \n")
  print(diff.df)
  
  # sd summary
  df.to.show = reshape(df.summary[c("GROUP", "MAF.bin", "sd")], direction = "wide", idvar = "GROUP", timevar = "MAF.bin")
  rownames(df.to.show) <- df.to.show$GROUP
  df.to.show = df.to.show[,2:dim(df.to.show)[2]]
  
  #TODO RAWNAK check how difference works with multiple groups
  diff.df = data.frame(diff(as.matrix(df.to.show)))
  
  colnames(df.to.show) = unlist(get.maf.bins.names())
  colnames(diff.df)    = unlist(get.maf.bins.names())
  printf("\n STANDART DEVIATION: \n")
  cat("=======================================================\n")
  print(df.to.show)
  
  printf("\n SD DIFFERENCE: (HA - controls): \n")
  print(diff.df)
}

eoh.outliers.analysis = function(in.df){
  # analysis of topright outliers 
  df.HA  = in.df[in.df$GROUP == "HA", 1:7]
  df.HA$EoH.controls = in.df$EoH[df$GROUP == "controls"]
  df.HA  = df.HA[df.HA$EoH >0.12,]
  df.HA  = df.HA[df.HA$EoH.controls >0.12,]
  
  # png("EoH_cc2.png")
  ggplot(df.HA, aes(EoH, EoH.controls, colour = MAF.bin)) + geom_point(shape = ".")
  ggsave("EoH_cc01.pdf")
  
  write.table(df.HA$SNP, file ="EoH_snp_012.txt", sep = "\t",
              row.names = FALSE, col.names = FALSE, quote=FALSE)
  
  EoH.snp.annotated = read.table("EoH_outliers_analysis/EoH_snp_annotation.txt", header =TRUE, as.is=T)
  gene.annotation = unique(EoH.snp.annotated[,c("Uploaded_variation", "Gene")])
  gene.annotation = gene.annotation[gene.annotation$Gene !="-",]
  colnames(gene.annotation) = c("SNP", "Gene")
  df.HA.annotated = base::merge(x = df.HA, y=gene.annotation, all.x = TRUE)
  df.HA.annotated = na.omit(df.HA.annotated)
  df.HA.annotated = unique(df.HA.annotated[with(df.HA.annotated, order(CHR)), c("CHR", "Gene")])
  write.table(df.HA.annotated, file ="outliers_gene_names.txt", sep = "\t",
                           row.names = FALSE, col.names = TRUE, quote=FALSE)
  #plot
  chr.position = data.frame(table(df.HA.annotated$CHR))
  colnames(chr.position) = c("CHR", "count")
  chr.position = chr.position[with(chr.position, order(-count)),]
  plot.count = chr.position[1:6,]
  levels(plot.count$CHR) = c(levels(plot.count$CHR),"other")
  plot.count = rbind(plot.count, c("other", sum(chr.position[7:nrow(chr.position), "count"])))
  plot.count[, "count"] = sapply(plot.count[, "count"], as.numeric)
  
  lbls <- paste("CHR ", plot.count$CHR, "\n", plot.count$count, sep="")
  png("pie.png")
  pie(plot.count$count, labels = lbls, radius = 1, col=rainbow(length(lbls)))
  dev.off() 
  
}
eoh.max.dif.analysis = function(df){
  df.plot = df[df$GROUP == "HA", c("CHR", "SNP", "EoH", "EoHOMOZ", "EoHOMOZ.major")]
  df.plot$EoH.controls = df$EoH[df$GROUP == "controls"]
  df.plot$EoHOMOZ.controls = df$EoHOMOZ[df$GROUP == "controls"]
  df.plot$EoHOMOZ.major.controls = df$EoHOMOZ.major[df$GROUP == "controls"]
  df.plot$EoH.dif = df.plot$EoH - df.plot$EoH.controls
  df.plot$EoHOMOZ.dif = df.plot$EoHOMOZ - df.plot$EoHOMOZ.controls
  df.plot$EoHOMOZ.major.dif = df.plot$EoHOMOZ.major - df.plot$EoHOMOZ.major.controls
  df.plot = base::merge(x = df.plot, y = plink.MAF, by = c("CHR", "SNP"))
  
  ggplot(df.plot, aes(MAF, EoH.dif)) + geom_point(shape = ".") + labs(title = "EoH: HA - controls")
  ggsave("EOH_dif_for_maf.png")
  ggplot(df.plot, aes(MAF, EoHOMOZ.dif)) + geom_point(shape = ".") + labs(title = "EoHOMOZ: HA - controls")
  ggsave("EoHOMOZ_dif_for_maf.png")
  ggplot(df.plot, aes(MAF, EoHOMOZ.major.dif)) + geom_point(shape = ".") + labs(title = "EoHOMOZ.major: HA - controls")
  ggsave("EoHOMOZ_major_dif_for_maf.png")
  
  #plot HLA SNPs
  snp.filter.path = "/projects/cgstudies/HA_GWAS_2017/Healthy_Aging_GWAS_2017_ABW-P01/Analysis/HAS_Het_Analysis_OLGA/Heterozygosity_tests/NWE_data/annotation/HLA.snp"
  snp.filter = read.table(sprintf("%s", snp.filter.path), header =TRUE, as.is=T)
  snp.filter$HLA = 1
  df.plot.cutoff = base::merge(x = df.plot.cutoff, y = snp.filter, by = c("CHR", "SNP"), all.x = TRUE)
  df.plot.cutoff$HLA[is.na(df.plot.cutoff$HLA)] = 0
  df.plot.cutoff$HLA=as.factor(df.plot.cutoff$HLA)
  
  ggplot(df.plot.cutoff, aes(MAF, EoH.dif, colour = HLA)) + geom_point(shape = ".") +
    labs(title = "EoH: HA - controls") + scale_color_manual(values=c("grey", "red"))
  ggsave("dif_for_maf_HLA.png")
  
  #plot SNP type
  snp.filter.path = "/projects/cgstudies/HA_GWAS_2017/Healthy_Aging_GWAS_2017_ABW-P01/Analysis/HAS_Het_Analysis_OLGA/Heterozygosity_tests/NWE_data/annotation/synomymous.snp"
  snp.filter = read.table(sprintf("%s", snp.filter.path), header =TRUE, as.is=T)
  snp.filter$SYN = 1
  df.plot.cutoff = base::merge(x = df.plot.cutoff, y = snp.filter, by = c("CHR", "SNP"), all.x = TRUE)
  df.plot.cutoff$SYN[is.na(df.plot.cutoff$SYN)] = 0
  df.plot.cutoff$SYN=as.factor(df.plot.cutoff$SYN)
  
  snp.filter.path = "/projects/cgstudies/HA_GWAS_2017/Healthy_Aging_GWAS_2017_ABW-P01/Analysis/HAS_Het_Analysis_OLGA/Heterozygosity_tests/NWE_data/annotation/nonsynonymous.snp"
  snp.filter = read.table(sprintf("%s", snp.filter.path), header =TRUE, as.is=T)
  snp.filter$NONSYN = 1
  df.plot.cutoff = base::merge(x = df.plot.cutoff, y = snp.filter, by = c("CHR", "SNP"), all.x = TRUE)
  df.plot.cutoff$NONSYN[is.na(df.plot.cutoff$NONSYN)] = 0
  df.plot.cutoff$NONSYN=as.factor(df.plot.cutoff$NONSYN)
  
  snp.filter.path = "/projects/cgstudies/HA_GWAS_2017/Healthy_Aging_GWAS_2017_ABW-P01/Analysis/HAS_Het_Analysis_OLGA/Heterozygosity_tests/NWE_data/annotation/intronic.snp"
  snp.filter = read.table(sprintf("%s", snp.filter.path), header =TRUE, as.is=T)
  snp.filter$INTR = 1
  df.plot.cutoff = base::merge(x = df.plot.cutoff, y = snp.filter, by = c("CHR", "SNP"), all.x = TRUE)
  df.plot.cutoff$INTR[is.na(df.plot.cutoff$INTR)] = 0
  df.plot.cutoff$INTR=as.factor(df.plot.cutoff$INTR)
  
  snp.filter.path = "/projects/cgstudies/HA_GWAS_2017/Healthy_Aging_GWAS_2017_ABW-P01/Analysis/HAS_Het_Analysis_OLGA/Heterozygosity_tests/NWE_data/annotation/UTR.snp"
  snp.filter = read.table(sprintf("%s", snp.filter.path), header =TRUE, as.is=T)
  snp.filter$UTR = 1
  df.plot.cutoff = base::merge(x = df.plot.cutoff, y = snp.filter, by = c("CHR", "SNP"), all.x = TRUE)
  df.plot.cutoff$UTR[is.na(df.plot.cutoff$UTR)] = 0
  df.plot.cutoff$UTR=as.factor(df.plot.cutoff$UTR)
  
  df.plot.cutoff$type = "other"
  df.plot.cutoff$type[df.plot.cutoff$INTR == 1] = "intron"
  df.plot.cutoff$type[df.plot.cutoff$UTR == 1] = "UTR"
  df.plot.cutoff$type[df.plot.cutoff$SYN == 1] = "SYN"
  df.plot.cutoff$type[df.plot.cutoff$NONSYN == 1] = "NONSYN"
  
  ggplot(df.plot.cutoff, aes(MAF, EoH.dif, colour = type)) + geom_point(shape = ".") +
    labs(title = "EoH: HA - controls") + scale_color_manual(values=c("lightyellow", "blue", "grey", "green","red"))
  ggsave("dif_for_maf_type.pdf")
}
eoh.regression.analysis = function(df){
  df.plot = df[, c("CHR", "SNP", "GROUP", "EoH")]
  
  snp.filter.path = "/projects/cgstudies/HA_GWAS_2017/Healthy_Aging_GWAS_2017_ABW-P01/Analysis/HAS_Het_Analysis_OLGA/Heterozygosity_tests/NWE_data/annotation/HLA.snp"
  snp.filter = read.table(sprintf("%s", snp.filter.path), header =TRUE, as.is=T)
  snp.filter$HLA = 1
  df.plot = base::merge(x = df.plot, y = snp.filter, by = c("CHR", "SNP"), all.x = TRUE)
  df.plot$HLA[is.na(df.plot$HLA)] = 0
  df.plot$HLA=as.factor(df.plot$HLA)
  
  snp.filter.path = "/projects/cgstudies/HA_GWAS_2017/Healthy_Aging_GWAS_2017_ABW-P01/Analysis/HAS_Het_Analysis_OLGA/Heterozygosity_tests/NWE_data/annotation/synomymous.snp"
  snp.filter = read.table(sprintf("%s", snp.filter.path), header =TRUE, as.is=T)
  snp.filter$SYN = 1
  df.plot = base::merge(x = df.plot, y = snp.filter, by = c("CHR", "SNP"), all.x = TRUE)
  df.plot$SYN[is.na(df.plot$SYN)] = 0
  df.plot$SYN=as.factor(df.plot$SYN)
  
  snp.filter.path = "/projects/cgstudies/HA_GWAS_2017/Healthy_Aging_GWAS_2017_ABW-P01/Analysis/HAS_Het_Analysis_OLGA/Heterozygosity_tests/NWE_data/annotation/nonsynonymous.snp"
  snp.filter = read.table(sprintf("%s", snp.filter.path), header =TRUE, as.is=T)
  snp.filter$NONSYN = 1
  df.plot = base::merge(x = df.plot, y = snp.filter, by = c("CHR", "SNP"), all.x = TRUE)
  df.plot$NONSYN[is.na(df.plot$NONSYN)] = 0
  df.plot$NONSYN=as.factor(df.plot$NONSYN)
  
  snp.filter.path = "/projects/cgstudies/HA_GWAS_2017/Healthy_Aging_GWAS_2017_ABW-P01/Analysis/HAS_Het_Analysis_OLGA/Heterozygosity_tests/NWE_data/annotation/intronic.snp"
  snp.filter = read.table(sprintf("%s", snp.filter.path), header =TRUE, as.is=T)
  snp.filter$INTR = 1
  df.plot = base::merge(x = df.plot, y = snp.filter, by = c("CHR", "SNP"), all.x = TRUE)
  df.plot$INTR[is.na(df.plot$INTR)] = 0
  df.plot$INTR=as.factor(df.plot$INTR)
  
  snp.filter.path = "/projects/cgstudies/HA_GWAS_2017/Healthy_Aging_GWAS_2017_ABW-P01/Analysis/HAS_Het_Analysis_OLGA/Heterozygosity_tests/NWE_data/annotation/UTR.snp"
  snp.filter = read.table(sprintf("%s", snp.filter.path), header =TRUE, as.is=T)
  snp.filter$UTR = 1
  df.plot = base::merge(x = df.plot, y = snp.filter, by = c("CHR", "SNP"), all.x = TRUE)
  df.plot$UTR[is.na(df.plot$UTR)] = 0
  df.plot$UTR=as.factor(df.plot$UTR)
  
  df.plot$type = "other"
  df.plot$type[df.plot$INTR == 1] = "intron"
  df.plot$type[df.plot$UTR == 1] = "UTR"
  df.plot$type[df.plot$SYN == 1] = "SYN"
  df.plot$type[df.plot$NONSYN == 1] = "NONSYN"
  
  save(df.plot, file = sprintf("%s/saved_data/%sdf.plot.RData", config$project.dir, config$subgroup.pref))
  }

library(lattice)
#https://github.com/etal/cnvkit/blob/master/cnvlib/scatter.py
#https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html for zscore
#https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R
manhattan.plot<-function(chr, pos, pvalue, 
                         sig.level=NA, annotate=NULL, ann.default=list(),
                         should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                         xlab="Chromosome", ylab=expression(-log[10](p-value)),
                         col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }
  
  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                      col=NULL, fontface=NULL, fontsize=NULL, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) +1 ], 
                              fill=lfills[(i-2) %% length(lfills) +1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  #reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      logp=round((pvalue),thin.logp.places), 
      pos=round(genpos,thin.pos.places), 
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- (pvalue)
  }
  rm(pos, pvalue)
  gc()
  
  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax+posmin)/2+posshift),
                 labels=levels(chr), 
                 ticks=F, rot=0,
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    }
    else {
      axis.default(side=side,...);
    }
  }
  
  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
    #A$ylim=c(0,maxy);
    A$ylim=c(min(y)+0.1,max(y)+0.1);
    A;
  }
  
  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings, 
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             #panel.abline(h=-log10(sig.level), lty=2);
             panel.abline(h=sig.level, lty=2);
             panel.abline(h=-sig.level, lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}
