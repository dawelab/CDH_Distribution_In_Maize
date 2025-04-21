library(ggplot2)
library(ggpubr)

#Load the data with all of the CDH calls
GROUPS <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumCompleteWholeGenomePCBChr.csv")
colnames(GROUPS)

SUB <- subset(GROUPS, Data_Source != "Romero-Navarro_etal_2017" & Maize_Type != "Inbred")

write.table(SUB$Name, "HighCoverageLines.txt", row.names = FALSE, quote = FALSE)

#Load the environmental data for all the GPS points from WorldClim2 and FAO
ENV <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")
colnames(ENV)

#This brings together the calls and the environmental data
DATA <- merge(GROUPS, ENV, by=c("Name", "Data_Source", "Accession", "Group", "Maize_Type", "Ab10_Status", "B_Chrom_Status", "Latitude", "Longitude", "Altitude", "DTA_BLUP", "DTS_BLUP", "PH_BLUP", "PresTiller_BLUP"))

colnames(DATA) <- gsub("wc2.1_30s_", "", colnames(DATA))

#This filters to only lines with environmental data
DATA_SUB <- subset(DATA, is.na(DATA$Latitude) == FALSE & DATA$KMeans_Ab10 != "Ambiguous")

a <- ggplot(DATA_SUB, aes(x=KMeans_Ab10, y=elev, color = KMeans_Ab10)) +
  geom_jitter(alpha=0.5) +
  geom_boxplot(alpha=0.5) +
  labs(x="Ab10 Status\nN.S.", y= "Elevation") +
  scale_color_manual(values=c("Ab10" = "chartreuse4", "N10" = "grey40")) +
  theme_classic() +
  theme(legend.position="none")
a

DATA_SUB <- subset(DATA, is.na(DATA$Latitude) == FALSE & DATA$KMeans_K10L2 != "Ambiguous")

b <- ggplot(DATA_SUB, aes(x=KMeans_K10L2, y=elev, color = KMeans_K10L2)) +
  geom_jitter(alpha=0.5) +
  geom_boxplot(alpha=0.5) +
  labs(x="K10L2 Status\np=0.03", y= "Elevation") +
  scale_color_manual(values=c("K10L2" = "red", "N10" = "grey40")) +
  theme_classic() +
  theme(legend.position="none")
b

DATA_SUB <- subset(DATA, is.na(DATA$Latitude) == FALSE & DATA$KMeans_BChrom != "Ambiguous")

DATA_SUB$KMeans_BChrom <- factor(DATA_SUB$KMeans_BChrom, levels=c("Yes", "No"))

c <- ggplot(DATA_SUB, aes(x=KMeans_BChrom, y=elev, color = KMeans_BChrom)) +
  geom_jitter(alpha=0.5) +
  geom_boxplot(alpha=0.5) +
  labs(x="B Chr. Status\nN.S.", y= "Elevation") +
  scale_color_manual(values=c("Yes" = "blue", "No" = "grey40")) +
  theme_classic() +
  theme(legend.position="none")
c


CopyNum <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/6_ClassifyAllData/BChr_ContExperimental_CopyNumber.csv")

#This brings together the calls and the environmental data
DATA2 <- merge(CopyNum, ENV, by=c("Name", "Data_Source", "Accession", "Group", "Maize_Type", "Ab10_Status", "B_Chrom_Status", "Latitude", "Longitude", "Altitude", "DTA_BLUP", "DTS_BLUP", "PH_BLUP", "PresTiller_BLUP"))

colnames(DATA2) <- gsub("wc2.1_30s_", "", colnames(DATA2))


pdf("Elevation_Scatter.pdf", height = 2, width=2.1)
d <- ggplot(DATA2, aes(x=Mean, y=elev)) +
  geom_jitter(alpha=0.5) +
  labs(x="B Chr.\nPseudo Copy Number\nN.S.", y= "Elevation") +
  #scale_color_manual(values=c("Yes" = "blue", "No" = "grey40")) +
  theme_classic() +
  xlim(0,0.21) +
  theme(legend.position="none")
d
dev.off()

pdf("Elevation_BoxPlots.pdf", height = 2, width=5)
ggarrange(a,b,c, nrow = 1)
dev.off()


