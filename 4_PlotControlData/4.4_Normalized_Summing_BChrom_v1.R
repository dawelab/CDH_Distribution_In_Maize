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

MERGE_BChrom_RPM <- vroom::vroom("Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.BChrom.RPM.txt")

#28 columns were classed as logicals because they have all NAs, these are the ones that had 0 coverage

BINS <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Bins_NoOverlap.table.xlsx")

#This drops any column that had an average coverage of 0, resulting in all values being NA
DT <- as.data.table(MERGE_BChrom_RPM)
DT <- DT[,which(unlist(lapply(DT, function(x)!all(is.na(x))))),with=F]
MERGE_BChrom_RPM <- as.data.frame(DT)

#This removes any tag with a MAPQ less than 20
MERGE_BChrom_RPM_FILT <- subset(MERGE_BChrom_RPM, MAPQ >= 20)

###This selects only the controls
#This selects only the control data
GROUPS_CONTROLS <- subset(GROUPS, GROUPS$Data_Source == "Dawe_Lab_1" | GROUPS$Data_Source == "Dawe_Lab_2")

#This drops blanks
GROUPS_CONTROLS <- subset(GROUPS_CONTROLS, Name != "BLANK.1.DC1" & Name != "BLANK.2.DC1" & Name != "BLANK.1.DC2" & Name != "NSL-2833_B-Chrom.2.DC2" & Name != "B542C_L289_B-Chrom.1.DC2")
nrow(GROUPS_CONTROLS)

#This identifies all the index numbers for names that appear in the GROUPS control only name fields
col.num <- which(colnames(MERGE_BChrom_RPM_FILT) %in% GROUPS_CONTROLS$Name)
#This appends the first 6 columns to that list
col.num <- append(col.num, values = c(1:6))
#This selects only the columns that appear in the above list
MERGE_BChrom_RPM_FILT_CONTROLS <- MERGE_BChrom_RPM_FILT[,sort(c(col.num))]
ncol(MERGE_BChrom_RPM_FILT_CONTROLS)

###################### This calculates total percent coverage on the Ab10 haplotype

#This calculates the sums of percents for each line
TotalSums_temp <- MERGE_BChrom_RPM_FILT_CONTROLS %>%
  summarise(across(c(7:ncol(MERGE_BChrom_RPM_FILT_CONTROLS)), ~ sum(., na.rm = TRUE)))

#This reshapes the new dataframe
#This corrects the column names
TotalSums_temp1 <- as.data.frame(t(TotalSums_temp))

#This makes the existing rownames (i.e. the line names) a column
TotalSums_temp1$Name <- row.names(TotalSums_temp1)

#This changes the column names
colnames(TotalSums_temp1) <- c("RPM", "Name")

#This merges the dataframe with Ab10 status 
GROUPS_CONTROLS_STAT <- GROUPS_CONTROLS[,c("Name", "B_Chrom_Status")]
TotalSums <- merge(TotalSums_temp1, GROUPS_CONTROLS_STAT, by = "Name")

ggplot(TotalSums, aes(x=B_Chrom_Status, y=RPM, fill=B_Chrom_Status, color=B_Chrom_Status)) +
  geom_violin() +
  geom_jitter(color = "black")+
  scale_color_viridis(option="viridis", discrete= TRUE) +
  scale_fill_viridis(option="viridis", discrete= TRUE) +
  ggtitle("RPM Tags Mapped to the BChromosome") +
  xlab("B Chromosome") +
  ylab("RPM") +
  theme(plot.title = element_text(hjust=0.5, size = 20), plot.subtitle = element_text(hjust=0.5, size = 15)) +
  guides(fill="none", color = "none")
ggsave("RPM_BChrom.png")

