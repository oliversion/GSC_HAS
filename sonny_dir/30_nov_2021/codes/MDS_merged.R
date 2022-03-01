setwd("/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/working")

options(bitmapType='cairo')
data<- read.table(file="MDS_merged1.mds",header=TRUE)
#data = data4
pop<- read.table(file="pop_file.txt",header=TRUE)
#datafile<- merge(data,pop,by=c("IID","FID"))

pop2 = pop[,c(2,3)]
datafile<- merge(data,pop2,by=c("IID"))


## Plot MDS 1 vs. MDS 2
pdf("MDS_1_2_plink.pdf",width=7,height=7)
colors <- c("black","green","red",470,"blue")
pchs = c(3,1,1,1,1)
reordered_groups <- factor(datafile$pop, levels = c("HAS","EUR","ASN","AMR","AFR"))
plot(datafile[,4], datafile[,5], col = colors[reordered_groups], pch=pchs[reordered_groups], 
     xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),
     main="MDS using PLINK: MDS1 vs. MDS2",
     xlab="MDS1",
     ylab="MDS2")
legend("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), pch=c(3,1,1,1,1), col = colors, cex=1)
dev.off()


## Plot MDS 1 vs. MDS 3
pdf("MDS_1_3_plink.pdf",width=7,height=7)
colors <- c("black","green","red",470,"blue")
pchs = c(3,1,1,1,1)
reordered_groups <- factor(datafile$pop, levels = c("HAS","EUR","ASN","AMR","AFR"))
plot(datafile[,4], datafile[,6], col = colors[reordered_groups], pch=pchs[reordered_groups], 
     xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),
     main="MDS using PLINK: MDS1 vs. MDS3",
     xlab="MDS1",
     ylab="MDS3")
legend("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), pch=c(3,1,1,1,1), col = colors, cex=1)
dev.off()


## Plot MDS 2 vs. MDS 3
pdf("MDS_2_3_plink.pdf",width=7,height=7)
colors <- c("black","green","red",470,"blue")
pchs = c(3,1,1,1,1)
reordered_groups <- factor(datafile$pop, levels = c("HAS","EUR","ASN","AMR","AFR"))
plot(datafile[,5], datafile[,6], col = colors[reordered_groups], pch=pchs[reordered_groups], 
     xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),
     main="MDS using PLINK: MDS2 vs. MDS3",
     xlab="MDS2",
     ylab="MDS3")
legend("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), pch=c(3,1,1,1,1), col = colors, cex=1)
dev.off()



