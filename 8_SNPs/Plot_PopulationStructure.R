library(ggplot2)

setwd("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/8_SNPs")

GROUPS <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/6.5_Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumComplete.csv")

#This loads all the eigen ventor values
Ab10_PCA <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10_PCA.eigenvec")
B_PCA <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChr_PCA.eigenvec")
K10L2_PCA <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/K10L2_PCA.eigenvec")

#This loads the percent of variation explained by each PC
Ab10_Val <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10_PCA.eigenval", header=FALSE)
Ab10_Val$Perc <- (Ab10_Val$V1/sum(Ab10_Val$V1))*100

K10L2_Val <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/K10L2_PCA.eigenval", header=FALSE)
K10L2_Val$Perc <- (K10L2_Val$V1/sum(K10L2_Val$V1))*100

B_Val <- read.delim("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChr_PCA.eigenval", header=FALSE)
B_Val$Perc <- (B_Val$V1/sum(B_Val$V1))*100

#This removes individuals without GPS coordinates so that the map and the PCA are comparable
GROUPS <- subset(GROUPS, is.na(Latitude)==FALSE)

############################################################################################################################################


#This merges the eigen vectors and the groups information 
MERGE_Ab10 <- merge(GROUPS, Ab10_PCA, by.x= "Name", by.y = "IID")

#This drops any samples without a PCA value likely because it was dropped in filtering
MERGE_Ab10 <- subset(MERGE_Ab10, is.na(PC1)==FALSE)

#This drops any ambigious samples and all control samples which have a Kmeans_Ab10 status of NA
MERGE_Ab10 <- subset(MERGE_Ab10, is.na(KMeans_Ab10)==FALSE, KMeans_Ab10 != "Ambiguous")

#This plots the points
pdf("Ab10_PCA.pdf", height=3, width=3)
ggplot() +
  geom_point(data=MERGE_Ab10[MERGE_Ab10$KMeans_Ab10 != "Ab10", ], aes(x=PC1, y=PC2, color=KMeans_Ab10), size=1) +
  geom_point(data=MERGE_Ab10[MERGE_Ab10$KMeans_Ab10 != "N10", ], aes(x=PC1, y=PC2, color=KMeans_Ab10), size=1) +
  scale_color_manual(values=c("Ab10"="chartreuse4", "N10"="grey")) +
  labs(x=paste0("PC1 (", round(Ab10_Val[[1,2]]), "%)"), y=paste0("PC2 (", round(Ab10_Val[[2,2]]), "%)")) +
  theme_classic() +
  xlim(-0.02,0.03) +
  theme(legend.position = "none")
dev.off()

#I am plotting this because Ab10 is associated with PC10
pdf("Ab10_PCA_2by10.pdf", height=3, width=3)
ggplot() +
  geom_point(data=MERGE_Ab10[MERGE_Ab10$KMeans_Ab10 != "Ab10", ], aes(x=PC2, y=PC10, color=KMeans_Ab10), size=1) +
  geom_point(data=MERGE_Ab10[MERGE_Ab10$KMeans_Ab10 != "N10", ], aes(x=PC2, y=PC10, color=KMeans_Ab10), size=1) +
  scale_color_manual(values=c("Ab10"="chartreuse4", "N10"="grey")) +
  labs(x=paste0("PC2 (", round(Ab10_Val[[2,2]]), "%)"), y=paste0("PC10 (", round(Ab10_Val[[10,2]]), "%)")) +
  theme_classic() +
  theme(legend.position = "none")
dev.off()


############################################################################################################################################

#This merges the eigen vectors and the groups information 
MERGE_K10L2 <- merge(GROUPS, K10L2_PCA, by.x= "Name", by.y = "IID")

#This drops any samples without a PCA value likely because it was dropped in filtering
MERGE_K10L2 <- subset(MERGE_K10L2, is.na(PC1)==FALSE)

#This drops any ambigious samples and all control samples which have a Kmeans_K10L2 status of NA
MERGE_K10L2 <- subset(MERGE_K10L2, is.na(KMeans_K10L2)==FALSE, KMeans_K10L2 != "Ambiguous")

#This plots the points
pdf("K10L2_PCA.pdf", height=3, width=3)
ggplot() +
  geom_point(data=MERGE_K10L2[MERGE_K10L2$KMeans_K10L2 != "K10L2", ], aes(x=PC1, y=PC2, color=KMeans_K10L2), size =1) +
  geom_point(data=MERGE_K10L2[MERGE_K10L2$KMeans_K10L2 != "N10", ], aes(x=PC1, y=PC2, color=KMeans_K10L2), size =1) +
  scale_color_manual(values=c("K10L2"="red", "N10"="grey")) +
  labs(x=paste0("PC1 (", round(K10L2_Val[[1,2]]), "%)"), y=paste0("PC2 (", round(K10L2_Val[[2,2]]), "%)")) +
  theme_classic() +
  xlim(-0.02,0.03) +
  theme(legend.position = "none")
dev.off()

############################################################################################################################################

#This merges the eigen vectors and the groups information 
MERGE_B <- merge(GROUPS, B_PCA, by.x= "Name", by.y = "IID")

#This drops any samples without a PCA value likely because it was dropped in filtering
MERGE_B <- subset(MERGE_B, is.na(PC1)==FALSE)

#This drops any ambigious samples and all control samples which have a Kmeans_B status of NA
MERGE_B <- subset(MERGE_B, is.na(KMeans_BChrom)==FALSE, KMeans_BChrom != "Ambiguous")

#This plots the points
pdf("BChr_PCA.pdf", height=3, width=3)
ggplot() +
  geom_point(data=MERGE_B[MERGE_B$KMeans_BChrom != "Yes", ], aes(x=PC1, y=PC2, color=KMeans_BChrom), size =1) +
  geom_point(data=MERGE_B[MERGE_B$KMeans_BChrom != "No", ], aes(x=PC1, y=PC2, color=KMeans_BChrom), size =1) +
  scale_color_manual(values=c("Yes"="blue", "No"="grey")) +
  labs(x=paste0("PC1 (", round(B_Val[[1,2]]), "%)"), y=paste0("PC2 (", round(B_Val[[2,2]]), "%)")) +
  theme_classic() +
  xlim(-0.02,0.03) +
  theme(legend.position = "none")
dev.off()
