setwd("/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/working/")
options(bitmapType='cairo')
data = read.table("HASpc.txt", header=T)
pop<- read.table(file="pop_file.txt",header=TRUE)
#datafile<- merge(data,pop,by=c("IID","FID"))

pop2 = pop[,c(2,3)]
datafile<- merge(data,pop2,by=c("IID"))

## Plot PCA 1 vs. PCA 2
pdf("MDS_1_2_king.pdf",width=7,height=7)
colors <- c("black","green","red",470,"blue")
pchs = c(3,1,1,1,1)
reordered_groups <- factor(datafile$pop, levels = c("HAS","EUR","ASN","AMR","AFR"))
plot(datafile[,7], datafile[,8], col = colors[reordered_groups], pch=pchs[reordered_groups], 
     xlim=c(-0.1,0.2),
     ylim=c(-0.15,0.1),
     main="MDS plot: Component 1 vs. 2",
     xlab="MDS 1",
     ylab="MDS 2")
legend("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), pch=c(3,1,1,1,1), col = colors, cex=1)
dev.off()