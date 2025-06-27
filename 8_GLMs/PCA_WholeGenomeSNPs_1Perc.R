#!/usr/bin/Rscript
library(stringr)
library(vcfR)
library(ggplot2)

setwd("/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/")
#This reads in the VCF file
VCF <- read.vcfR("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/CoordinateOnly_Ab10BChrom.chr1to9.filt4.1Perc.vcf.gz")
#This extracts only the genotypes
print("#####Extracting genotypes")
GT<- extract.gt(VCF)

#This converts the genotypes to an additive coding so that it can be used in a PCA
print("#####Converting genotypes to additive coding")
convert <- function(x) {
  y <- str_split((x), ":")[[1]][1]
  y <- sub("0/1", "1", y)
  y <- sub("1/0", "1", y)
  y <- sub("1/1", "2", y)
  y <- sub("0/0", "0", y)
  x <<-y
}
GT <- apply(GT, c(1, 2), convert)

print("#####Transposing data")
for_PCA <- t(GT)

#This creates the distance matrix and then performs the PCA
print("#####Creating distance matrix")
dist <- dist(for_PCA)
print("#####Performing PCA")
pca <- cmdscale(dist, k=10, eig = TRUE)

#This calculates the percent of the variation they account for 
print("#####Making scree plot")
TOT <- sum(pca$eig)
pca$perc_var <- (pca$eig/TOT)*100 

#This plots the eigen vector values. This tells you how much variability is accounted for by each PCA
pdf("PCA_ScreePlot_1Perc.pdf", height = 5, width = 5)
plot(pca$perc_var)
dev.off()

#This converts the PCA to a data frame
print("#####Converting PCA output to a data frame")
df_out <- as.data.frame(matrix(nrow = length(rownames(for_PCA)), ncol = 11))
names(df_out) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Name")
df_out$PC1 <- pca$points[,1]
df_out$PC2 <- pca$points[,2]
df_out$PC3 <- pca$points[,3]
df_out$PC4 <- pca$points[,4]
df_out$PC5 <- pca$points[,5]
df_out$PC6 <- pca$points[,6]
df_out$PC7 <- pca$points[,7]
df_out$PC8 <- pca$points[,8]
df_out$PC9 <- pca$points[,9]
df_out$PC10 <- pca$points[,10]
df_out$Name <-  c(as.character(row.names(for_PCA)))

#Write out the PCA results
print("#####Writing the PCA dataframe to a file")
write.csv(df_out, "WholeGenomeSNP_PCA_1Perc.csv", row.names=FALSE, quote=FALSE)

#This plots it
print("#####Plotting PCA")
pdf("PCA1_v_PCA2_1Perc.pdf", height = 8, width = 15)
ggplot(df_out, aes(x=PC1, y=PC2)) +
  geom_point(alpha = 0.8, size = 2 ) +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()

pdf("PCA1_v_PCA3_1Perc.pdf", height = 8, width = 15)
ggplot(df_out, aes(x=PC1, y=PC3)) +
  geom_point(alpha = 0.8, size = 2 ) +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()

pdf("PCA1_v_PCA3_1Perc.pdf", height = 8, width = 15)
ggplot(df_out, aes(x=PC1, y=PC4)) +
  geom_point(alpha = 0.8, size = 2 ) +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()

pdf("PCA2_v_PCA3_1Perc.pdf", height = 8, width = 15)
ggplot(df_out, aes(x=PC2, y=PC3)) +
  geom_point(alpha = 0.8, size = 2 ) +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()