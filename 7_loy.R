# import libraries====================================================================
library(tictoc)
library(ggplot2)

#====================================================================
#      READ DATA
# ====================================================================
report.dir = ""
sample.map = read.csv(sprintf("%s/%s", report.dir, "Sample_Map.txt"), header =T, as.is=T, sep = "\t")
# raw.Y.region = read.csv(sprintf("%s/%s", report.dir, "out_Y_region.txt"), header =F, as.is=T, sep = "\t")
raw.Y.region = read.csv(sprintf("%s/%s", report.dir, "Out.txt"), header =F, as.is=T, sep = "\t")
colnames(raw.Y.region) = c("SNP.Name",	"Sample.ID",	"Allele1",	"Allele2",	"GC.Score",	"Sample.Name",	
                           "Sample.Group",	"SNP.Index",	"Chr",	"Position", "SNP",	"Theta",	
                           "R",	"X",	"Y", "X.Raw",	"Y.Raw",	"BAF",	"LLR")
unique(raw.Y.region$Sample.Name) # 74 samples
df.plot = ddply(raw.Y.region, .(as.factor(raw.Y.region$Sample.Name)), nrow) 
unique(df.plot$V1) # all samples have 948 SNPs at male specific Y Chr
# temp = raw.Y.region[raw.Y.region$Sample.Name== "C365_ENS",]
samles.males = sample.map$Name[sample.map$Gender=="Male"]

raw.Y.region.males = raw.Y.region[raw.Y.region$Sample.Name %in% samles.males,]
unique(raw.Y.region.males$Sample.Name)

# aggregate(LLR ~ Sample.Name, data = raw.Y.region.males, summary)
summary.LLR = raw.Y.region.males%>%
  group_by(Sample.Name)%>% 
  summarise(Mean=mean(LLR, na.rm = T), Max=max(LLR, na.rm = T), Min=min(LLR, na.rm = T),
            Median=median(LLR, na.rm = T), Std=sd(LLR, na.rm = T))

#====================================================================
#      PLOTS
# ====================================================================
#example B Allele Freq	Log R Ratio
sample.Y = raw.Y.region.males[raw.Y.region.males$Sample.Name== "C278_DES",] # C278_DES normal; LOY HA16-692
p1 = ggplot(sample.Y, aes(x=LLR)) +
  geom_histogram(colour="black", fill="white", bins = 60) + 
  xlab("Log R Ratio")+
  ggtitle("No evidence for mosaicism: C278_DES")
p2 = ggplot(sample.Y, aes(x=BAF)) +
  geom_histogram(colour="black", fill="white", bins = 60) +
  xlab("B Allele Frequency") +
  ggtitle("No evidence for mosaicism: C278_DES")
sample.Y = raw.Y.region.males[raw.Y.region.males$Sample.Name== "HA16-692",]
p3 = ggplot(sample.Y, aes(x=LLR)) +
    geom_histogram(colour="black", fill="white", bins = 60) + 
    xlab("Log R Ratio")+
    ggtitle("LOY -0.8: HA16-692")
p4 = ggplot(sample.Y, aes(x=BAF)) +
    geom_histogram(colour="black", fill="white", bins = 60) +
  xlab("B Allele Frequency")+
    ggtitle("LOY -0.8: HA16-692")
file.name = sprintf("~/output/12_LOY/%s.%s", "LOY_example", "pdf")
ggsave(filename=file.name, plot=multiplot(p1, p2, cols=2), height = 8, width = 14)
  

