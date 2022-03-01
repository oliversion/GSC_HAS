setwd("/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/working")
options(bitmapType='cairo', rgl.printRglwidget=TRUE)
data<- read.table(file="PCA_merged1.eigenvec",header=T, comment.char = "")

pop<- read.table(file="pop_file.txt",header=TRUE)
#datafile<- merge(data,pop,by=c("IID","FID"))

pop2 = pop[,c(2,3)]
datafile<- merge(data,pop2,by=c("IID"))


## Plot PCA 1 vs. PCA 2
pdf("PCA_1_2_plink.pdf",width=7,height=7)

colors <- c("black","green","red",470,"blue")
pchs = c(3,1,1,1,1)
reordered_groups <- factor(datafile$pop, levels = c("HAS","EUR","ASN","AMR","AFR"))
plot(datafile[,3], datafile[,4], col = colors[reordered_groups], pch=pchs[reordered_groups], 
     xlim=c(-0.1,0.2),
     ylim=c(-0.15,0.1),
     main="PCA using PLINK: PC1 vs. PC2",
     xlab="PC1",
     ylab="PC2")
legend("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), pch=c(3,1,1,1,1), col = colors, cex=1)

dev.off()


## Plot PCA 1 vs. PCA 3
pdf("PCA_1_3_plink.pdf",width=7,height=7)

colors <- c("black","green","red",470,"blue")
pchs = c(3,1,1,1,1)
reordered_groups <- factor(datafile$pop, levels = c("HAS","EUR","ASN","AMR","AFR"))
plot(datafile[,3], datafile[,5], col = colors[reordered_groups], pch=pchs[reordered_groups], 
     xlim=c(-0.1,0.2),
     ylim=c(-0.15,0.1),
     main="PCA using PLINK: PC1 vs. PC3",
     xlab="PC1",
     ylab="PC3")
legend("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), pch=c(3,1,1,1,1), col = colors, cex=1)

dev.off()


## Plot PCA 2 vs. PCA 3
pdf("PCA_2_3_plink.pdf",width=7,height=7)

colors <- c("black","green","red",470,"blue")
pchs = c(3,1,1,1,1)
reordered_groups <- factor(datafile$pop, levels = c("HAS","EUR","ASN","AMR","AFR"))
plot(datafile[,4], datafile[,5], col = colors[reordered_groups], pch=pchs[reordered_groups], 
     xlim=c(-0.1,0.2),
     ylim=c(-0.15,0.1),
     main="PCA using PLINK: PC2 vs. PC3",
     xlab="PC2",
     ylab="PC3")
legend("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), pch=c(3,1,1,1,1), col = colors, cex=1)

dev.off()


## plot 3D scatterplot

# library(rgl)
# plot3d( 
#   x=datafile[,3], y=datafile[,4], z=datafile[,5], 
#   col = colors[reordered_groups], 
#   type = 's', 
#   radius = .002,
#   xlab="PCA1",
#   ylab="PCA2",
#   zlab="PCA3"
#   )
# legend3d("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), col = colors, pch=16, cex=1, inset=c(0.02))

