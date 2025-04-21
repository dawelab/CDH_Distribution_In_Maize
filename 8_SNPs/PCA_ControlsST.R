library(stringr)

#I manually removed the header from the file output by Filter_Higgins_SNPs.R

VCF <- read.table("/Volumes/Transcend/AllData_v5.Ab10Shared.ControlsST.filt.nohead.vcf", header=TRUE)

convert <- function(x) {
  y <- str_split((x), ":")[[1]][1]
  y <- sub("0/1", "1", y)
  y <- sub("1/0", "1", y)
  y <- sub("1/1", "2", y)
  y <- sub("0/0", "0", y)
  if(grepl(".", y, fixed = TRUE)==TRUE) {
    y<- 0
  }
  x <<-y
}

VCF[,-c(1:9)] <- apply(VCF[,-c(1:9)], c(1, 2), convert)

#This reformats the data for compatibility with StAMPP
VCF <- VCF[,-c(1:9)]


#transpose data 
for_PCA <- t(VCF)

#This creates the distance matrix and then performs the PCA
dist <- dist(for_PCA)
pca <- cmdscale(dist, k=10, eig = TRUE)

#This calculates the percent of the variation they account for 
TOT <- sum(pca$eig)
pca$perc_var <- (pca$eig/TOT)*100 

#This plots the eigen vector values. This tells you how much variability is accounted for by each PCA
pdf("PCA_ScreePlot.pdf", height = 5, width = 5)
plot(pca$perc_var)
dev.off()

#This converts the PCA to a data frame
df_out <- as.data.frame(matrix(nrow = length(rownames(for_PCA)), ncol = 9))
names(df_out) <- c("PC1", "PC2", "PC3", "PC4", "Name")
df_out$PC1 <- pca$points[,1]
df_out$PC2 <- pca$points[,2]
df_out$PC3 <- pca$points[,3]
df_out$PC4 <- pca$points[,4]
df_out$Name <-  c(as.character(row.names(for_PCA)))
df_out$Chr10Type <- c("K10L2", "K10L2", "K10L2", "K10L2", "K10L2", "K10L2", "K10L2", "Ab10-III", "Ab10-III", "Ab10-III", "Ab10-III", "Ab10-III", "N10", "N10", "Ab10-II", "Ab10-II", "Ab10-I", "Ab10-I", "Ab10-I", "Ab10-I", "N10", "Ab10-II", "Ab10-II", "Ab10-II", "K10L2", "Ab10-II", "Ab10-II", "K10L2", "Ab10-II", "Ab10-I", "Ab10-II", "Ab10-II", "Ab10-III", "Ab10-III", "Ab10-III", "Ab10-III", "Ab10-III", "N10", "N10", "N10", "Ab10-III", "Ab10-III", "Ab10-III", "Ab10-III", "Ab10-III", "Ab10-II", "Ab10-II", "Ab10-II", "Ab10-II", "Ab10-II", "N10", "N10", "N10", "Ab10-I", "Ab10-I", "Ab10-I", "Ab10-I", "Ab10-I", "N10", "N10", "N10")

#This plots it
pdf("PCA1_v_PCA2.pdf", height = 8, width = 15)
ggplot(df_out, aes(x=PC1, y=PC2, color=Chr10Type)) +
  geom_point(alpha = 0.8, size = 2 ) +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()

pdf("PCA1_v_PCA3.pdf", height = 8, width = 15)
ggplot(df_out, aes(x=PC1, y=PC3, color=Chr10Type)) +
  geom_point(alpha = 0.8, size = 2 ) +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()

pdf("PCA1_v_PCA3.pdf", height = 8, width = 15)
ggplot(df_out, aes(x=PC1, y=PC4, color=Chr10Type)) +
  geom_point(alpha = 0.8, size = 2 ) +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()

pdf("PCA2_v_PCA3.pdf", height = 8, width = 15)
ggplot(df_out, aes(x=PC2, y=PC3, color=Chr10Type)) +
  geom_point(alpha = 0.8, size = 2 ) +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()
