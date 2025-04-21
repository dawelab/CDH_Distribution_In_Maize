GROUPS <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/6.5_Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumComplete.csv")

K10L2 <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/6_ClassifyAllData/K10L2_ContExperimental_CopyNumber.csv")
K10L2 <- K10L2[,c("Name", "Mean")]
colnames(K10L2) <- c("Name", "K10L2_Pseudo_Copy_Number")
Ab10 <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/6_ClassifyAllData/Ab10_ContExperimental_CopyNumber.csv")
Ab10 <- Ab10[,c("Name", "Mean")]
colnames(Ab10) <- c("Name", "Ab10_Pseudo_Copy_Number")
BChrom <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/6_ClassifyAllData/BChr_ContExperimental_CopyNumber.csv")
BChrom <- BChrom[,c("Name", "Mean")]
colnames(BChrom) <- c("Name", "BChrom_Pseudo_Copy_Number")

GROUPS_NEW <- merge(GROUPS, Ab10, all=TRUE)
GROUPS_NEW <- merge(GROUPS_NEW, K10L2, all=TRUE)
GROUPS_NEW <- merge(GROUPS_NEW, BChrom, , all=TRUE)

write.csv(GROUPS_NEW, "/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/6.5_Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumComplete_v2.csv", row.names = FALSE, quote=FALSE) 
