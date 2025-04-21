library(ggplot2)

K10L2_BCOMB <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/K10L2_Model/K10L2_TagSumPlusTagDensAllSamples.csv")
SCG_BCOMB <- vroom::vroom("/Volumes/Transcend/SCG_TagSumPlusTagDensAllSamples.csv")

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

MISCLASS <-  c("W23_AB10-I.11.DC1", "W23_AB10-I.13.DC1", "W23_AB10-II.36.DC1", "W23_N10.14.DC1")

#This drops misclassified controls
COV_CONT <- COV_CONT[-c(which(COV_CONT$Name %in% MISCLASS)),]

#This merges with the groups file
COV_CONT_GROUPS <- merge(COV_CONT, GROUPS, by ="Name")


COV_CONT_GROUPS_K10L2 <- subset(COV_CONT_GROUPS, Ab10_Status == "K10L2" | Ab10_Status == "N10")

ggplot(data=COV_CONT_GROUPS_K10L2, aes(x=Ab10_Status, y=Mean)) +
  geom_jitter(height = 0) 
ggsave("K10L2_CopyNumber_Controls.png")

COV_EXP <- COV[-c(grep("DC", COV$Name)),]
COV_EXP_GROUPS <- merge(COV_EXP, GROUPS, by = "Name")

COV_EXP_GROUPS_K10L2 <- subset(COV_EXP_GROUPS, KMeans_K10L2 == "K10L2")

COV_ALL_GROUPS <- rbind(COV_CONT_GROUPS_K10L2, COV_EXP_GROUPS_K10L2)
COV_ALL_GROUPS$Ab10_Status <- factor(COV_ALL_GROUPS$Ab10_Status, levels = c("N10", "K10L2", "Unknown"))

ggplot(data=COV_ALL_GROUPS, aes(x=Ab10_Status, y=Mean)) +
  geom_jitter(height=0)
ggsave("K10L2_CopyNumber_Experimental.png")

ggplot(data=COV_ALL_GROUPS, aes(x=Ab10_Status, y=Mean)) +
  geom_jitter(height=0) +
  ylim(0,0.2)
ggsave("K10L2_CopyNumber_Experimental_axisfix.png")


