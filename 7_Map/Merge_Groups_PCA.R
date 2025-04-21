
GROUPS <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/6_ClassifyAllData/6.5_Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumComplete.csv")

PCA <- read.table("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10_PCA.eigenvec", header=TRUE)

PCA <- PCA[,-c(1)]
colnames(PCA) <- c("Name", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

BOTH <- merge(GROUPS, PCA, by="Name")

write.csv(BOTH, "~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumCompleteWholeGenomePCAb10.csv", row.names = FALSE, quote = FALSE)

##################################################################

PCA <- read.table("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/K10L2_PCA.eigenvec", header=TRUE)

PCA <- PCA[,-c(1)]
colnames(PCA) <- c("Name", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

BOTH <- merge(GROUPS, PCA, by="Name")

write.csv(BOTH, "~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumCompleteWholeGenomePCK10L2.csv", row.names = FALSE, quote = FALSE)

##################################################################

PCA <- read.table("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChr_PCA.eigenvec", header=TRUE)

PCA <- PCA[,-c(1)]
colnames(PCA) <- c("Name", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

BOTH <- merge(GROUPS, PCA, by="Name")

write.csv(BOTH, "~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumCompleteWholeGenomePCBChr.csv", row.names = FALSE, quote = FALSE)

