library(ggplot2)


Ab10_BCOMB <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/Ab10Model/Ab10_TagSumPlusTagDensAllSamples.csv")

SCG_BCOMB <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/SCG_TagSumPlusTagDensAllSamples.csv")

GROUPS <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromComplete.csv")

Mean <- apply(SCG_BCOMB, 2, mean)
Median <- apply(SCG_BCOMB, 2, median)

SCG_Sum <- data.frame(Name=colnames(SCG_BCOMB), Mean= Mean, Median = Median)

######This one does the best 
COV <- data.frame(Name = c(NA), Mean = c(NA))
COV <- COV[-c(1),]

i=1
for(i in 1:ncol(Ab10_BCOMB)) {
  NAME <- colnames(Ab10_BCOMB)[i]
  EDIT <- Ab10_BCOMB[,which(colnames(Ab10_BCOMB) %in% NAME)] / SCG_Sum[which(SCG_Sum$Name %in% NAME),2]
  MEAN <- mean(EDIT[-c(1:2, 13:17, 19:27, 50),])
  LINE <- data.frame(Name=NAME, Mean = MEAN)
  COV <<- rbind(COV, LINE)
}

COV_CONT <- COV[grep("DC", COV$Name),]

#This drops misclassified samples and the Ab10 Unknown sample, and the B73 samples because they are misleading
MISCLASS <-  c("W23_AB10-I.11.DC1", "W23_AB10-I.13.DC1", "W23_AB10-II.36.DC1", "W23_N10.14.DC1", "PI-483314_K10L2.3.DC2", "B73_N10.1.DC2", "B73_N10.2.DC2", "B73_N10.3.DC2")

#This drops misclassified controls
COV_CONT <- COV_CONT[-c(which(COV_CONT$Name %in% MISCLASS)),]

#This defines the homozygous lines. 
# Removed "W23_AB10-II.22.DC1", because it was clearly not Hom
#HOM <- c("FFMM_Ab10-I_hom.1.DC2", "FFMM_Ab10-I_hom.2.DC2", "FFMM_Ab10-I_hom.3.DC2", "W23_AB10-I.1.DC1", "W23_AB10-I.36.DC1", "W23_AB10-I.17.DC1", "W23_AB10-I.38.DC1", "W23_AB10-I.5.DC1", "W23_AB10-I.8.DC1", "W23_AB10-I.3.DC1", "W23_AB10-I.4.DC1", "W23_AB10-II.34.DC1", "W23_AB10-II.29.DC1", "W23_AB10-II.14.DC1", "W23_AB10-II.6.DC1", "W23_AB10-II.28.DC1", "W23_AB10-II.33.DC1", "W23_AB10-II.1.DC1", "W23_AB10-II.11.DC1", "W23_AB10-II.17.DC1", "W23_AB10-II.25.DC1", "W23_AB10-II.21.DC1", "W23_AB10-II.24.DC1", "W23_AB10-II.19.DC1", "W23_AB10-II.10.DC1", "W23_AB10-II.5.DC1", "W23_AB10-II.30.DC1")

#This drops samples with an unknown coy number
Drop <- grep(".DC1", COV_CONT$Name)
COV_CONT <- COV_CONT[-c(Drop),]

#This is a more stringent assesment of Hom
HOM <- c("FFMM_Ab10-I_hom.1.DC2", "FFMM_Ab10-I_hom.2.DC2", "FFMM_Ab10-I_hom.3.DC2")

COV_CONT_HOM <- COV_CONT[which(COV_CONT$Name %in% HOM),]
COV_CONT_HOM$Copy <- "2"

COV_CONT_HET <- COV_CONT[-c(which(COV_CONT$Name %in% HOM)),]

GROUPS_N10 <- subset(GROUPS, Ab10_Status == "N10" | Ab10_Status == "K10L2")
COV_CONT_HET <- COV_CONT_HET[-c(which(COV_CONT_HET$Name %in% GROUPS_N10$Name)),]
COV_CONT_HET$Copy <- "1"

GROUPS_TrueN10 <- subset(GROUPS, Ab10_Status == "N10")
COV_CONT_N10 <- COV_CONT[c(which(COV_CONT$Name %in% GROUPS_TrueN10$Name)),]
COV_CONT_N10$Copy <- "0"

COV_CONT <- rbind(COV_CONT_HOM, COV_CONT_HET, COV_CONT_N10)

COV_CONT_GROUPS <- merge(COV_CONT, GROUPS, by ="Name")

#This appends all of the experimental samples 
COV_Groups <- merge(COV, GROUPS, by ="Name")

#This selects only B positive samples that are not controls
COV_GroupsAb10 <- subset(COV_Groups, KMeans_Ab10 == "Ab10" & Data_Source != "Dawe_Lab_1" & Data_Source != "Dawe_Lab_2")

COV_GroupsAb10$Copy <- "Experimental"

ALL <- rbind(COV_CONT_GROUPS, COV_GroupsAb10)
ALL$CDH <- "Ab10"

write.csv(ALL, "Ab10_ContExperimental_CopyNumber.csv", row.names = FALSE, quote = FALSE)
