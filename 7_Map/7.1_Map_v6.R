install.packages("maps")
install.packages("mapdata")
library(maps)
library(mapdata)
library(ggplot2)

GROUPS <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/7.5_Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumComplete_Edit.csv", header = TRUE)

########################################################
#Ab10
########################################################

Ab10_Yes <- subset(GROUPS, GROUPS$KMeans_Ab10 == "Ab10")
Ab10_No <- subset(GROUPS, GROUPS$KMeans_Ab10 == "N10")

III <-  subset(GROUPS, GROUPS$RF_Ab10Type == "Ab10-III")
II <-  subset(GROUPS, GROUPS$RF_Ab10Type == "Ab10-II")
I <-  subset(GROUPS, GROUPS$RF_Ab10Type == "Ab10-I")
Ambig <- subset(GROUPS, GROUPS$RF_Ab10Type == "Ambiguous")

png("Ab10_Map.png")
map("world", xlim = c(-140, -30), ylim = c(-60, 50))
points(Ab10_No$Longitude, Ab10_No$Latitude, pch=20, col="darkgrey", cex=0.2)
points(Ab10_Yes$Longitude, Ab10_Yes$Latitude, pch=20, col="blue", cex=0.5)
dev.off()

png("Ab10-I_TypeMap_Map.png")
map("world", xlim = c(-140, -30), ylim = c(-60, 50))
points(I$Longitude, I$Latitude, pch=20, col="darkorange", cex=0.5)
dev.off()

png("Ab10-II_TypeMap_Map.png")
map("world", xlim = c(-140, -30), ylim = c(-60, 50))
points(II$Longitude, II$Latitude, pch=20, col="darkgreen", cex=1)
dev.off()

png("Ab10-III_TypeMap_Map.png")
map("world", xlim = c(-140, -30), ylim = c(-60, 50))
points(III$Longitude, III$Latitude, pch=20, col="purple", cex=1)
dev.off()

png("Ab10-Ambig_TypeMap_Map.png")
map("world", xlim = c(-140, -30), ylim = c(-60, 50))
points(Ambig$Longitude, Ambig$Latitude, pch=20, col="black", cex=0.5)
dev.off()

########################################################
#K10L2
########################################################

K10L2_Yes <- subset(GROUPS, GROUPS$KMeans_K10L2 == "K10L2")
K10L2_No <- subset(GROUPS, GROUPS$KMeans_K10L2 == "N10")

png("K10L2_Map.png")
map("world", xlim = c(-140, -30), ylim = c(-60, 50))
points(K10L2_No$Longitude, K10L2_No$Latitude, pch=20, col="darkgrey", cex=0.2)
points(K10L2_Yes$Longitude, K10L2_Yes$Latitude, pch=20, col="red", cex=0.5)
dev.off()

########################################################
#K10L2
########################################################

BChrom_Yes <- subset(GROUPS, GROUPS$KMeans_BChrom == "Yes")
BChrom_No <- subset(GROUPS, GROUPS$KMeans_BChrom == "No")

png("BChrom_Map.png")
map("world", xlim = c(-140, -30), ylim = c(-60, 50))
points(BChrom_No$Longitude, BChrom_No$Latitude, pch=20, col="darkgrey", cex=0.2)
points(BChrom_Yes$Longitude, BChrom_Yes$Latitude, pch=20, col="turquoise4", cex=0.5)
dev.off()

B_Chrom_1 <- subset(GROUPS, GROUPS$KMeans_BChrom == "Yes" & BChrom_Bin == 1)
B_Chrom_2 <- subset(GROUPS, GROUPS$KMeans_BChrom == "Yes" & BChrom_Bin == 2)
B_Chrom_3 <- subset(GROUPS, GROUPS$KMeans_BChrom == "Yes" & BChrom_Bin == 3)
B_Chrom_4 <- subset(GROUPS, GROUPS$KMeans_BChrom == "Yes" & BChrom_Bin == 4)
B_Chrom_5 <- subset(GROUPS, GROUPS$KMeans_BChrom == "Yes" & BChrom_Bin == 5)
B_Chrom_6 <- subset(GROUPS, GROUPS$KMeans_BChrom == "Yes" & BChrom_Bin == 6)
B_Chrom_7 <- subset(GROUPS, GROUPS$KMeans_BChrom == "Yes" & BChrom_Bin == 7)
B_Chrom_8 <- subset(GROUPS, GROUPS$KMeans_BChrom == "Yes" & BChrom_Bin == 8)

png("BChrom_Map.png")
map("world", xlim = c(-140, -30), ylim = c(-60, 50))
points(BChrom_No$Longitude, BChrom_No$Latitude, pch=20, col="darkgrey", cex=0.2)
points(B_Chrom_1$Longitude, B_Chrom_1$Latitude, pch=20, col="#D5EECD", cex=0.5)
points(B_Chrom_2$Longitude, B_Chrom_2$Latitude, pch=20, col="#BAE4BD", cex=0.5)
points(B_Chrom_3$Longitude, B_Chrom_3$Latitude, pch=20, col="#9CD8B8", cex=0.5)
points(B_Chrom_4$Longitude, B_Chrom_4$Latitude, pch=20, col="#7BCCC4", cex=0.5)
points(B_Chrom_5$Longitude, B_Chrom_5$Latitude, pch=20, col="#59B9CF", cex=0.5)
points(B_Chrom_6$Longitude, B_Chrom_6$Latitude, pch=20, col="#3C9FC8", cex=0.5)
points(B_Chrom_7$Longitude, B_Chrom_7$Latitude, pch=20, col= "#227FB6", cex=0.5)
points(B_Chrom_8$Longitude, B_Chrom_8$Latitude, pch=20, col="#08589E", cex=0.5)
dev.off()
