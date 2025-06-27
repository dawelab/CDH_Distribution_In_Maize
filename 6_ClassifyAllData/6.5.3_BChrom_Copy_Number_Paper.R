library(ggplot2)

BChrom_BCOMB <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/BChrom/BChrom_TagSumPlusTagDensAllSamples.csv")
SCG_BCOMB <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/SCG_TagSumPlusTagDensAllSamples.csv")

GROUPS <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromComplete.csv")

Mean <- apply(SCG_BCOMB, 2, mean)
Median <- apply(SCG_BCOMB, 2, median)

SCG_Sum <- data.frame(Name=colnames(SCG_BCOMB), Mean= Mean, Median = Median)

######This one does the best 
COV <- data.frame(Name = c(NA), Mean = c(NA))
COV <- COV[-c(1),]

i=1
for(i in 1:ncol(BChrom_BCOMB)) {
  NAME <- colnames(BChrom_BCOMB)[i]
  EDIT <- BChrom_BCOMB[,which(colnames(BChrom_BCOMB) %in% NAME)] / SCG_Sum[which(SCG_Sum$Name %in% NAME),2]
  MEAN <- mean(EDIT[-c(4, 107, 2, 5, 87),])
  LINE <- data.frame(Name=NAME, Mean = MEAN)
  COV <<- rbind(COV, LINE)
}

COV_CONT <- COV[grep("DC", COV$Name),]

MISCLASS <-  c("NSL-2833_B-Chrom.2.DC2", "B542C_L289_B-Chrom.1.DC2")

#This drops misclassified controls
COV_CONT <- COV_CONT[-c(which(COV_CONT$Name %in% MISCLASS)),]

#This merges with the groups file
COV_CONT_GROUPS <- merge(COV_CONT, GROUPS, by ="Name")
COV_CONT_GROUPS_NoB <- subset(COV_CONT_GROUPS, B_Chrom_Status == "No")
COV_CONT_GROUPS_NoB$Copy <- "0"


COV_CONT_GROUPS_B <- subset(COV_CONT_GROUPS, B_Chrom_Status == "Yes")
Low <- c("NSL-2833_B-Chrom.3.DC2", "B542C_L289_B-Chrom.2.DC2", "B542C_L289_B-Chrom.3.DC2", "B542C_L289_B-Chrom.4.DC2")
COV_CONT_GROUPS_BL <- COV_CONT_GROUPS_B[which(COV_CONT_GROUPS_B$Name %in% Low),]
COV_CONT_GROUPS_BL$Copy <- "Unknown"
COV_CONT_GROUPS_BH <- COV_CONT_GROUPS_B[-c(which(COV_CONT_GROUPS_B$Name %in% Low)),]
COV_CONT_GROUPS_BH$Copy <- "Unknown"

COV_CONT_GROUPS_B <- rbind(COV_CONT_GROUPS_NoB, COV_CONT_GROUPS_BL, COV_CONT_GROUPS_BH)

#This appends all of the experimental samples 
COV_Groups <- merge(COV, GROUPS, by ="Name")

#This selects only B positive samples that are not controls
COV_GroupsB <- subset(COV_Groups, KMeans_BChrom == "Yes" & Data_Source != "Dawe_Lab_1" & Data_Source != "Dawe_Lab_2")

COV_GroupsB$Copy <- "Experimental"

ALL <- rbind(COV_CONT_GROUPS_B, COV_GroupsB)
ALL$CDH <- "BChr"

write.csv(ALL, "BChr_ContExperimental_CopyNumber.csv", row.names = FALSE, quote = FALSE)
