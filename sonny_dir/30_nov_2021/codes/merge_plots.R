setwd("/projects/cgstudies/HAS_Diversity_OLGA/data/sonny_dir/30-nov-2021/1__/")
king_pca_raw = read.table("HAS_pca.txt", header=T)
plink_pca_raw = read.table("PCA_merged1.eigenvec", header=T, comment.char = "")

pop<- read.table(file="pop_file.txt",header=TRUE)
pop2 = pop[,c(2,3)]
king_pca<- merge(king_pca_raw,pop2,by=c("IID"))
plink_pca = merge(plink_pca_raw, pop2, by=c("IID"))


colors <- c("black","green","red",470,"blue")
pchs = c(3,1,1,1,1)
reordered_groups_king <- factor(king_pca$pop, levels = c("HAS","EUR","ASN","AMR","AFR"))
reordered_groups_plink <- factor(plink_pca$pop, levels = c("HAS","EUR","ASN","AMR","AFR"))

jpeg("PCA_merged.jpg", height = 720)
par(mfrow = c(3,2))
?jpeg

### KING PC1 vs. PC2
plot(king_pca[,7], king_pca[,8], col = colors[reordered_groups_king], pch=pchs[reordered_groups_king], 
     xlim=c(-0.075,0.10),
     ylim=c(-0.1,0.095),
     main="PCA using KING: PC1 vs. PC2",
     xlab="PC1",
     ylab="PC2"
     )
#abline(v=-0.005, col="gray", lty=2) #h=c(-0.01,0.005), 
legend("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), pch=c(3,1,1,1,1), col = colors, cex=1)


### PLINK PC1 vs. PC2
plot(plink_pca[,3], plink_pca[,4], col = colors[reordered_groups_plink], pch=pchs[reordered_groups_plink], 
     xlim=c(-0.075,0.10),
     ylim=c(-0.1,0.095),
     main="PCA using PLINK: PC1 vs. PC2",
     xlab="PC1",
     ylab="PC2"
     )
legend("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), pch=c(3,1,1,1,1), col = colors, cex=1)


### KING PC1 vs. PC3
plot(king_pca[,7], king_pca[,9], col = colors[reordered_groups_king], pch=pchs[reordered_groups_king], 
     xlim=c(-0.075,0.10),
     ylim=c(-0.1,0.095),
     main="PCA using KING: PC1 vs. PC3",
     xlab="PC1",
     ylab="PC3")
legend("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), pch=c(3,1,1,1,1), col = colors, cex=1)


### PLINK PC1 vs. PC3
plot(plink_pca[,3], plink_pca[,5], col = colors[reordered_groups_plink], pch=pchs[reordered_groups_plink], 
     xlim=c(-0.075,0.10),
     ylim=c(-0.1,0.095),
     main="PCA using PLINK: PC1 vs. PC3",
     xlab="PC1",
     ylab="PC3")
legend("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), pch=c(3,1,1,1,1), col = colors, cex=1)


### KING PC2 vs. PC3
plot(king_pca[,8], king_pca[,9], col = colors[reordered_groups_king], pch=pchs[reordered_groups_king], 
     xlim=c(-0.075,0.10),
     ylim=c(-0.1,0.095),
     main="PCA using KING: PC2 vs. PC3",
     xlab="PC2",
     ylab="PC3")
legend("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), pch=c(3,1,1,1,1), col = colors, cex=1)




### PLINK PC2 vs. PC3
plot(plink_pca[,4], plink_pca[,5], col = colors[reordered_groups_plink], pch=pchs[reordered_groups_plink], 
     xlim=c(-0.075,0.10),
     ylim=c(-0.1,0.095),
     main="PCA using PLINK: PC2 vs. PC3",
     xlab="PC2",
     ylab="PC3")
legend("topright", legend = c("HAS","EUR","ASN","AMR","AFR"), pch=c(3,1,1,1,1), col = colors, cex=1)

dev.off()