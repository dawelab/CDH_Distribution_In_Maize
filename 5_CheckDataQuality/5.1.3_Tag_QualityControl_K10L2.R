library(ggplot2)
library(vroom)
library(tidyverse)
library(data.table)
library(readxl)
library(stringr)

########################################This loads and preps the data 
#This file comes from 1.7
GROUPS <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This file comes from 3.4
MERGE_K10L2_RPM <- vroom::vroom("Tassel_TagTaxaDist_AllData_v7_v_K10L2.K10L2.RPM.txt")
#8 columns were classed as logicals because they have all NAs, these are the ones that had 0 coverage

#This counts the number of non 0 tags in each line before MAPQ filtering
NUM_TAG <- apply(MERGE_K10L2_RPM[,-c(1:6)], 2, function(x) sum( x != 0))

#This removes any tag with a MAPQ less than 20
MERGE_K10L2_RPM_FILT <- subset(MERGE_K10L2_RPM, MAPQ >= 20)

#This counts the number of non 0 tags in each line after MAPQ filtering
NUM_TAG_FILT <- apply(MERGE_K10L2_RPM_FILT[,-c(1:6)], 2, function(x) sum( x != 0))

#This calculates the number of tags lost to filtering for each line 
QC <- as.data.frame(cbind(NUM_TAG, NUM_TAG_FILT))
#This calculates the number of tags that were removed due to low MAPQ
QC$Perc <- 100-(QC$NUM_TAG_FILT/QC$NUM_TAG)*100

QC$Name <- rownames(QC)
QC$Simp_Name <- NA
QC$Extra_Name <- NA
QC$OG_Name <- NA

#This section isolates the Romero Navarro Data and merges it with the groups.
GROUPS_RN <- subset(GROUPS, Data_Source == "Romero-Navarro_etal_2017")
QC_RN_temp <- QC[(grep("RN", QC$Name)),]

#This function splits the Romero Navarro names so they are compatible with the groups file
x <- 1
SPLIT <- function(x){
  QC_RN_temp[x,5] <<- paste(strsplit(QC_RN_temp[x,4], "[.]")[[1]][1], "RN", sep = ".")
  QC_RN_temp[x,6] <<- paste((strsplit(QC_RN_temp[x,4], "[.]")[[1]][2]), (strsplit(QC_RN_temp[x,4], "[.]")[[1]][3]), sep=".")
  QC_RN_temp[x,7] <<- QC_RN_temp[x,4]
}

lapply(c(1:nrow(QC_RN_temp)), SPLIT)

QC_RN <- merge(QC_RN_temp, GROUPS_RN, by.x ="Simp_Name", by.y = "Name")

#############Merging all of these separatly makes trouble shooting easier. 
#This section isolates the Swarts Data and merges it with the groups.
GROUPS_SW <- subset(GROUPS, Data_Source == "Swarts_etal_2017")
QC_SW <-merge(QC, GROUPS_SW, by = "Name")

#This section isolates the Romay Data and merges it with the groups.
GROUPS_RY <- subset(GROUPS, Data_Source == "Romay_etal_2013")
QC_RY <- merge(QC, GROUPS_RY, by = "Name")

#This section isolates the Romay Data and merges it with the groups.
GROUPS_DC <- subset(GROUPS, Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2")
QC_DC <- merge(QC, GROUPS_DC, by = "Name")

#This brings together all of the merged data
QC_v2 <- rbind(QC_DC, QC_RY, QC_SW, QC_RN)

#This plots the distribution of the percentage of tags lost to the MAPQ filtering by dataset
png(filename = "K10L2_PercentTagsFilteredForMAPQByData.png", height = 700, width = 500)
ggplot(data=QC_v2) +
  geom_histogram(aes(x=Perc, fill=Data_Source), bins = 100) +
  facet_wrap(~Data_Source, scales="free_y", ncol = 1) +
  xlab("Percent of Tags Filtered Out For MAPQ") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text = element_text(size = 15), legend.position = "none")
dev.off()

QC_v2_CONT <- subset(QC_v2, QC_v2$Data_Source == "Dawe_Lab_1" | QC_v2$Data_Source == "Dawe_Lab_2")

png(filename = "K10L2_PercentTagsFilteredForMAPQByK10L2Status.png", height = 1100, width = 500)
ggplot(data=QC_v2_CONT) +
  geom_histogram(aes(x=Perc, fill=Ab10_Status), bins = 100) +
  facet_wrap(~Ab10_Status, scales="fixed", ncol = 1) +
  xlim(0,100) +
  xlab("Percent of Tags Filtered Out For MAPQ") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text = element_text(size = 15), legend.position = "none")
dev.off()
