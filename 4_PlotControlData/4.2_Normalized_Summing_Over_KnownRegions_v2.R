library(ggplot2)
library(vroom)
library(tidyverse)
library(data.table)
library(stringr)
library(MASS)
library(reshape2)
library(viridis)
library(readxl)


#This sets the working directory
setwd("/Volumes/Transcend")

########################################This loads and preps the data 
GROUPS <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

MERGE_Ab10Hap_RPM <- vroom::vroom("Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.Ab10.RPM.txt")

#28 columns were classed as logicals because they have all NAs, these are the ones that had 0 coverage

BINS <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Bins_NoOverlap.table.xlsx")

#This drops any column that had an average coverage of 0, resulting in all values being NA
DT <- as.data.table(MERGE_Ab10Hap_RPM)
DT <- DT[,which(unlist(lapply(DT, function(x)!all(is.na(x))))),with=F]
MERGE_Ab10Hap_RPM <- as.data.frame(DT)

#This removes any tag with a MAPQ less than 20
MERGE_Ab10Hap_RPM_FILT <- subset(MERGE_Ab10Hap_RPM, MAPQ >= 20)

###This selects only the controls
#This selects only the control data
GROUPS_CONTROLS <- subset(GROUPS, GROUPS$Data_Source == "Dawe_Lab_1" | GROUPS$Data_Source == "Dawe_Lab_2")
GROUPS_CONTROLS$Name[grep("BLANK", GROUPS_CONTROLS$Name)]

#This drops lines that appear to be misclassified and or are unknown
GROUPS_CONTROLS <- subset(GROUPS_CONTROLS, Name != "W23_AB10-I.11.DC1" & Name != "W23_AB10-II.36.DC1" & Name != "W23_AB10-I.13.DC1" & Name != "W23_AB10-I.30.DC1" & Name != "W23_N10.14.DC1" & Name != "PI-483314_K10L2.3.DC2" & Name != "BLANK.1.DC1" & Name != "BLANK.2.DC1" & Name != "BLANK.1.DC2")
nrow(GROUPS_CONTROLS)

#This identifies all the index numbers for names that appear in the GROUPS control only name fields
col.num <- which(colnames(MERGE_Ab10Hap_RPM_FILT) %in% GROUPS_CONTROLS$Name)
#This appends the first 6 columns to that list
col.num <- append(col.num, values = c(1:6))
#This selects only the columns that appear in the above list
MERGE_Ab10Hap_RPM_FILT_CONTROLS <- MERGE_Ab10Hap_RPM_FILT[,sort(c(col.num))]
ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS)

#This assigns a region to each tag 
MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG <- MERGE_Ab10Hap_RPM_FILT_CONTROLS
MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG$Region <- "NA"

for(i in 1:nrow(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)) {
  VALUE <- MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,3]
  MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)] <- ifelse(VALUE >= 142472000 & VALUE <= 146699300, 'TR1', MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)])
  MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)] <- ifelse(VALUE >= 148964528 & VALUE <= 149082763, 'trkin', MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)])
  MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)] <- ifelse(VALUE >= 150656000 & VALUE <= 153145000, 'TR1', MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)])
  MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)] <- ifelse(VALUE >= 152050000 & VALUE <= 156350000, 'Shared region', MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)])
  MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)] <- ifelse(VALUE >= 158250000 & VALUE <= 166820000, 'Shared region', MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)])
  MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)] <- ifelse(VALUE >= 157485200 & VALUE <= 159356550, 'TR1', MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)])
  MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)] <- ifelse(VALUE >= 174433450 & VALUE <= 182846100, 'knob 180', MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)])
  MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)] <- ifelse(VALUE >= 189326066 & VALUE <= 190330226, 'kindr', MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)])
  MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)] <- ifelse(VALUE >= 190307761 & VALUE <= 191107857, 'kin10-like', MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG)])
}

summary(as.factor(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG$Region))

#This sums across the regions
RegionSums_temp <- MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG %>%
  group_by(Region)  %>%
  summarise(across(c(7:189), ~ sum(., na.rm = TRUE)))

#This reshapes the new dataframe
#This corrects the column names
RegionSums_temp1 <- as.data.frame(t(RegionSums_temp))
colnames(RegionSums_temp1) <- c("NA", "Shared Region", "TR1", "kin10-like", "kindr", "knob180", "trkin")
RegionSums_temp1 <- RegionSums_temp1[-c(1),]


#This makes the existing rownames (i.e. the line names) a column
RegionSums_temp1$Name <- row.names(RegionSums_temp1)

RegionSums_temp2 <- melt(RegionSums_temp1, id="Name")
colnames(RegionSums_temp2) <- c("Name", "Feature", "RPM")

#This merges the dataframe with Ab10 status 
GROUPS_CONTROLS_STAT <- GROUPS_CONTROLS[,c("Name", "Ab10_Status")]
RegionSums <- merge(RegionSums_temp2, GROUPS_CONTROLS_STAT, by = "Name")

#This removes all temporary variables
rm(RegionSums_temp)
rm(RegionSums_temp1)
rm(RegionSums_temp2)

#This converts the RPM column to a numeric. I am not sure why it was stored as a character
RegionSums$RPM <- as.numeric(RegionSums$RPM)

################## This plots all regions

#This sets the x axis order
level_order_all <- c("TR1", "Shared Region", "trkin", "knob180", "kindr", "kin10-like", "NA") 

ggplot(RegionSums, aes(x=factor(Feature, level = level_order_all), y=RPM, color=Ab10_Status, fill=Ab10_Status)) +
  geom_violin() +
  ggtitle("RPM Within Regions Of Known Importance") +
  xlab("Region of Known Importance") +
  ylab("RPM") +
  scale_color_viridis(option="viridis", discrete= TRUE) +
  scale_fill_viridis(option="viridis", discrete= TRUE) +
  geom_vline(xintercept = 1.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 2.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 3.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 4.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 5.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 6.5, linetype="dotted", color="darkgrey") +
  theme(plot.title = element_text(hjust=0.5, size = 20), plot.subtitle = element_text(hjust=0.5, size = 15)) +
  guides(fill=guide_legend(title="Chr10 Haplotype"), color= "none")

ggsave("RPM_Regions_Of_Known_Importance_All.png")


################## This plots only genes to make it easier to look at
RegionSumsGenes <- subset(RegionSums, Feature == "trkin" | Feature == "kindr" | Feature == "kin10-like")

level_order_genes <- c("trkin", "kindr", "kin10-like") 

ggplot(RegionSumsGenes, aes(x=factor(Feature, level = level_order_genes), y=RPM, color=Ab10_Status, fill=Ab10_Status)) +
  geom_violin() +
  ggtitle("RPM Within Regions Of Known Importance") +
  xlab("Region of Known Importance") +
  ylab("RPM") +
  scale_color_viridis(option="viridis", discrete= TRUE) +
  scale_fill_viridis(option="viridis", discrete= TRUE) +
  geom_vline(xintercept = 1.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 2.5, linetype="dotted", color="darkgrey") +
  theme(plot.title = element_text(hjust=0.5, size = 20), plot.subtitle = element_text(hjust=0.5, size = 15)) +
  guides(fill=guide_legend(title="Chr10 Haplotype"), color = "none")
ggsave("RPM_Regions_Of_Known_Importance_Genes.png")

###################### This calculates total percent coverage on the Ab10 haplotype

#This calculates the sums of percents for each line
TotalSums_temp <- MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG %>%
  summarise(across(c(7:189), ~ sum(., na.rm = TRUE)))

#This reshapes the new dataframe
#This corrects the column names
TotalSums_temp1 <- as.data.frame(t(TotalSums_temp))

#This makes the existing rownames (i.e. the line names) a column
TotalSums_temp1$Name <- row.names(TotalSums_temp1)

#This changes the column names
colnames(TotalSums_temp1) <- c("RPM", "Name")

#This merges the dataframe with Ab10 status 
GROUPS_CONTROLS_STAT <- GROUPS_CONTROLS[,c("Name", "Ab10_Status")]
TotalSums <- merge(TotalSums_temp1, GROUPS_CONTROLS_STAT, by = "Name")

ggplot(TotalSums, aes(x=Ab10_Status, y=RPM, fill=Ab10_Status, color=Ab10_Status)) +
  geom_violin() +
  scale_color_viridis(option="viridis", discrete= TRUE) +
  scale_fill_viridis(option="viridis", discrete= TRUE) +
  ggtitle("RPM Tags Mapped to the Ab10 Haplotype") +
  xlab("Chr10 Haplotype") +
  ylab("RPM") +
  theme(plot.title = element_text(hjust=0.5, size = 20), plot.subtitle = element_text(hjust=0.5, size = 15)) +
  guides(fill="none", color = "none")
ggsave("RPM_Ab10.png")


################## This plots everything EXCEPT the shared region summed together

#This calculates the sums of percents for each line
UniqueSums_temp <- subset(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG, Region != "Shared region")

#This checks that the shared region has been removed
summary(as.factor(UniqueSums_temp$Region))

#This sums each individual line
UniqueSums_temp1 <- UniqueSums_temp %>%
  summarise(across(c(7:189), ~ sum(., na.rm = TRUE)))

#This reshapes the new dataframe
#This corrects the column names
UniqueSums_temp2 <- as.data.frame(t(UniqueSums_temp1))

#This makes the existing rownames (i.e. the line names) a column
UniqueSums_temp2$Name <- row.names(UniqueSums_temp2)

#This changes the column names
colnames(UniqueSums_temp2) <- c("RPM", "Name")

#This merges the dataframe with Ab10 status 
GROUPS_CONTROLS_STAT <- GROUPS_CONTROLS[,c("Name", "Ab10_Status")]
UniqueSums <- merge(UniqueSums_temp2, GROUPS_CONTROLS_STAT, by = "Name")

ggplot(UniqueSums, aes(x=Ab10_Status, y=RPM, fill=Ab10_Status, color=Ab10_Status)) +
  geom_violin() +
  scale_color_viridis(option="viridis", discrete= TRUE) +
  scale_fill_viridis(option="viridis", discrete= TRUE) +
  ggtitle("RPM Tags Mapped to the Ab10 Specific Region") +
  xlab("Chr10 Haplotype") +
  ylab("RPM") +
  theme(plot.title = element_text(hjust=0.5, size = 20), plot.subtitle = element_text(hjust=0.5, size = 15)) +
  guides(fill="none", color = "none")
ggsave("RPM_Ab10Specific.png")



##########################################################################
########### This section does all the same plots without K10L2############
##########################################################################

RegionSumsSub <- subset(RegionSums, Ab10_Status != "K10L2")

#This sets the x axis order
level_order_all <- c("TR1", "Shared Region", "trkin", "knob180", "kindr", "kin10-like", "NA") 

ggplot(RegionSumsSub, aes(x=factor(Feature, level = level_order_all), y=RPM, color=Ab10_Status, fill=Ab10_Status)) +
  geom_violin() +
  ggtitle("RPM Within Regions Of Known Importance") +
  xlab("Region of Known Importance") +
  ylab("RPM") +
  scale_color_viridis(option="viridis", discrete= TRUE) +
  scale_fill_viridis(option="viridis", discrete= TRUE) +
  geom_vline(xintercept = 1.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 2.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 3.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 4.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 5.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 6.5, linetype="dotted", color="darkgrey") +
  theme(plot.title = element_text(hjust=0.5, size = 20), plot.subtitle = element_text(hjust=0.5, size = 15)) +
  guides(fill=guide_legend(title="Chr10 Haplotype"), color= "none")

ggsave("RPM_Regions_Of_Known_Importance_All_NoK10L2.png")


################## This plots only genes to make it easier to look at
RegionSumsSubGenes <- subset(RegionSumsSub, Feature == "trkin" | Feature == "kindr" | Feature == "kin10-like")

level_order_genes <- c("trkin", "kindr", "kin10-like") 

ggplot(RegionSumsSubGenes, aes(x=factor(Feature, level = level_order_genes), y=RPM, color=Ab10_Status, fill=Ab10_Status)) +
  geom_violin() +
  ggtitle("RPM Within Regions Of Known Importance") +
  xlab("Region of Known Importance") +
  ylab("RPM") +
  scale_color_viridis(option="viridis", discrete= TRUE) +
  scale_fill_viridis(option="viridis", discrete= TRUE) +
  geom_vline(xintercept = 1.5, linetype="dotted", color="darkgrey") +
  geom_vline(xintercept = 2.5, linetype="dotted", color="darkgrey") +
  theme(plot.title = element_text(hjust=0.5, size = 20), plot.subtitle = element_text(hjust=0.5, size = 15)) +
  guides(fill=guide_legend(title="Chr10 Haplotype"), color = "none")
ggsave("RPM_Regions_Of_Known_Importance_Genes_NoK10L2.png")

###################### This calculates total percent coverage on the Ab10 haplotype

#This calculates the sums of percents for each line
TotalSums_temp <- MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG %>%
  summarise(across(c(7:189), ~ sum(., na.rm = TRUE)))

#This reshapes the new dataframe
#This corrects the column names
TotalSums_temp1 <- as.data.frame(t(TotalSums_temp))

#This makes the existing rownames (i.e. the line names) a column
TotalSums_temp1$Name <- row.names(TotalSums_temp1)

#This changes the column names
colnames(TotalSums_temp1) <- c("RPM", "Name")

#This merges the dataframe with Ab10 status 
GROUPS_CONTROLS_STAT <- GROUPS_CONTROLS[,c("Name", "Ab10_Status")]
GROUPS_CONTROLS_STAT_SUB <- subset(GROUPS_CONTROLS_STAT, Ab10_Status != "K10L2")
TotalSums <- merge(TotalSums_temp1, GROUPS_CONTROLS_STAT_SUB, by = "Name")

ggplot(TotalSums, aes(x=Ab10_Status, y=RPM, fill=Ab10_Status, color=Ab10_Status)) +
  geom_violin() +
  scale_color_viridis(option="viridis", discrete= TRUE) +
  scale_fill_viridis(option="viridis", discrete= TRUE) +
  ggtitle("RPM Tags Mapped to the Ab10 Haplotype") +
  xlab("Chr10 Haplotype") +
  ylab("RPM") +
  theme(plot.title = element_text(hjust=0.5, size = 20), plot.subtitle = element_text(hjust=0.5, size = 15)) +
  guides(fill="none", color = "none")
ggsave("RPM_Ab10_NoK10L2.png")


################## This plots everything EXCEPT the shared region summed together

#This calculates the sums of percents for each line
UniqueSums_temp <- subset(MERGE_Ab10Hap_RPM_FILT_CONTROLS_REG, Region != "Shared region")

#This checks that the shared region has been removed
summary(as.factor(UniqueSums_temp$Region))

#This sums each individual line
UniqueSums_temp1 <- UniqueSums_temp %>%
  summarise(across(c(7:189), ~ sum(., na.rm = TRUE)))

#This reshapes the new dataframe
#This corrects the column names
UniqueSums_temp2 <- as.data.frame(t(UniqueSums_temp1))

#This makes the existing rownames (i.e. the line names) a column
UniqueSums_temp2$Name <- row.names(UniqueSums_temp2)

#This changes the column names
colnames(UniqueSums_temp2) <- c("RPM", "Name")

#This merges the dataframe with Ab10 status 
GROUPS_CONTROLS_STAT <- GROUPS_CONTROLS[,c("Name", "Ab10_Status")]
GROUPS_CONTROLS_STAT_SUB <- subset(GROUPS_CONTROLS_STAT, Ab10_Status != "K10L2")
UniqueSums <- merge(UniqueSums_temp2, GROUPS_CONTROLS_STAT_SUB, by = "Name")

ggplot(UniqueSums, aes(x=Ab10_Status, y=RPM, fill=Ab10_Status, color=Ab10_Status)) +
  geom_violin() +
  scale_color_viridis(option="viridis", discrete= TRUE) +
  scale_fill_viridis(option="viridis", discrete= TRUE) +
  ggtitle("RPM Tags Mapped to the Ab10 Specific Region") +
  xlab("Chr10 Haplotype") +
  ylab("RPM") +
  theme(plot.title = element_text(hjust=0.5, size = 20), plot.subtitle = element_text(hjust=0.5, size = 15)) +
  guides(fill="none", color = "none")
ggsave("RPM_Ab10Specific_N0K10L2.png")


