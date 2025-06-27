library(ggplot2)

BChrom_BCOMB <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/BChrom/BChrom_TagSumPlusTagDensAllSamples.csv")
SCG_BCOMB <- vroom::vroom("/Volumes/Transcend/SCG_TagSumPlusTagDensAllSamples.csv")

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
COV_CONT_GROUPS_BL$Copy <- "Low"
COV_CONT_GROUPS_BH <- COV_CONT_GROUPS_B[-c(which(COV_CONT_GROUPS_B$Name %in% Low)),]
COV_CONT_GROUPS_BH$Copy <- "High"

COV_CONT_GROUPS_B <- rbind(COV_CONT_GROUPS_NoB, COV_CONT_GROUPS_BL, COV_CONT_GROUPS_BH)

ggplot(data=COV_CONT_GROUPS_B, aes(x=B_Chrom_Status, y=Mean, color = Copy)) +
  geom_jitter(height = 0) + 
  geom_hline(yintercept = 0.02, color="blue", linetype='dotted') +
  geom_hline(yintercept = 0.04, color="blue", linetype='dotted') +
  geom_hline(yintercept = 0.06, color="blue", linetype='dotted') +
  geom_hline(yintercept = 0.08, color="blue", linetype='dotted') +
  geom_hline(yintercept = 0.1, color="blue", linetype='dotted') +
  geom_hline(yintercept = 0.12, color="blue", linetype='dotted')
ggsave("BChrom_CopyNumber_Controls.png")


COV_EXP <- COV[-c(grep("DC", COV$Name)),]
COV_EXP_GROUPS <- merge(COV_EXP, GROUPS, by = "Name")
COV_EXP_GROUPS$Copy <- "Unknown"

COV_EXP_GROUPS_B <- subset(COV_EXP_GROUPS, KMeans_BChrom == "Yes")

COV_ALL_GROUPS <- rbind(COV_CONT_GROUPS_B, COV_EXP_GROUPS_B)
COV_ALL_GROUPS$B_Chrom_Status <- factor(COV_ALL_GROUPS$B_Chrom_Status, levels = c("No", "Yes", "Unknown"))

ggplot(data=COV_ALL_GROUPS, aes(x=B_Chrom_Status, y=Mean, color = Model)) +
  geom_jitter(height=0)
ggsave("BChrom_CopyNumber_Experimental.png")

#This adds the final pseudo copy number to the Groups file
COV_EDIT <- COV
colnames(COV_EDIT) <- c("Name", "B_PseudoCopyNum")
COV_GROUPS <- merge(GROUPS, COV, by = "Name")

write.csv(COV_GROUPS, file="/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumComplete.csv", row.names = FALSE, quote = FALSE)

COV_ALL_GROUPS_EXP <- subset(COV_ALL_GROUPS, Data_Source != "Dawe_Lab_1" & Data_Source != "Dawe_Lab_2")

ggplot(data=COV_ALL_GROUPS_EXP, aes(x=Maize_Type, y=Mean, color = Model)) +
  geom_jitter(height=0, width = 0.25, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("BChrom_CopyNumber_Experimental_MaizeType.png")

ggplot(data=COV_ALL_GROUPS_EXP, aes(x=Latitude, y=Mean, color = Maize_Type)) +
  geom_point(alpha=0.5) +
  ylab("Pseudo Copy Number") 
ggsave("BChrom_CopyNumber_by_Latitude.png")

ggplot(data=COV_ALL_GROUPS_EXP, aes(x=Longitude, y=Mean, color = Maize_Type)) +
  geom_point(alpha=0.5) +
  ylab("Pseudo Copy Number")
ggsave("BChrom_CopyNumber_by_Longitude.png")


ggplot(data=COV_ALL_GROUPS_EXP, aes(x=Altitude, y=Mean, color = Maize_Type)) +
  
  geom_point(alpha=0.5) +
  ylab("Pseudo Copy Number") +
  xlim(0,4000) +
ggsave("BChrom_CopyNumber_by_Elevation.png")
