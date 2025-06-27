library(ggplot2)

K10L2_BCOMB <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/K10L2_Model/K10L2_TagSumPlusTagDensAllSamples.csv")

SCG_BCOMB <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/SCG_TagSumPlusTagDensAllSamples.csv")

GROUPS <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromComplete.csv")

Mean <- apply(SCG_BCOMB, 2, mean)
Median <- apply(SCG_BCOMB, 2, median)

SCG_Sum <- data.frame(Name=colnames(SCG_BCOMB), Mean= Mean, Median = Median)

######This one does the best 
COV <- data.frame(Name = c(NA), Mean = c(NA))
COV <- COV[-c(1),]

i=1
for(i in 1:ncol(K10L2_BCOMB)) {
  NAME <- colnames(K10L2_BCOMB)[i]
  EDIT <- K10L2_BCOMB[,which(colnames(K10L2_BCOMB) %in% NAME)] / SCG_Sum[which(SCG_Sum$Name %in% NAME),2]
  #This intentionally pulls these rows, which is different than the other models which are exclusing them
  MEAN <- mean(EDIT[c(7:10),])
  LINE <- data.frame(Name=NAME, Mean = MEAN)
  COV <<- rbind(COV, LINE)
}

COV_CONT <- COV[grep("DC", COV$Name),]

MISCLASS <-  c("W23_AB10-I.11.DC1", "W23_AB10-I.13.DC1", "W23_AB10-II.36.DC1", "W23_N10.14.DC1", "PI-483314_K10L2.3.DC2")

#This drops misclassified controls
COV_CONT <- COV_CONT[-c(which(COV_CONT$Name %in% MISCLASS)),]

#This is a more stringent assesment of Hom
HOM <- c("PI-483314_K10L2.3.DC2", "B73_K10L2.10.DC2", "B73_K10L2.4.DC2")

COV_CONT_HOM <- COV_CONT[which(COV_CONT$Name %in% HOM),]
COV_CONT_HOM$Copy <- "2"

COV_CONT_HET <- COV_CONT[-c(which(COV_CONT$Name %in% HOM)),]

GROUPS_K10L2 <- subset(GROUPS, Ab10_Status == "K10L2")
COV_CONT_HET <- COV_CONT_HET[c(which(COV_CONT_HET$Name %in% GROUPS_K10L2$Name)),]
COV_CONT_HET$Copy <- "1"

GROUPS_TrueN10 <- subset(GROUPS, Ab10_Status == "N10")
COV_CONT_N10 <- COV_CONT[c(which(COV_CONT$Name %in% GROUPS_TrueN10$Name)),]
COV_CONT_N10$Copy <- "0"

COV_CONT <- rbind(COV_CONT_HOM, COV_CONT_HET, COV_CONT_N10)

COV_CONT_GROUPS <- merge(COV_CONT, GROUPS, by ="Name")

#This appends all of the experimental samples 
COV_Groups <- merge(COV, GROUPS, by ="Name")

#This selects only B positive samples that are not controls
COV_GroupsK10L2 <- subset(COV_Groups, KMeans_K10L2 == "K10L2" & Data_Source != "Dawe_Lab_1" & Data_Source != "Dawe_Lab_2")

COV_GroupsK10L2$Copy <- "Experimental"

ALL <- rbind(COV_CONT_GROUPS, COV_GroupsK10L2)

ALL$CDH <- "K10L2"

write.csv(ALL, "K10L2_ContExperimental_CopyNumber.csv", row.names = FALSE, quote = FALSE)
