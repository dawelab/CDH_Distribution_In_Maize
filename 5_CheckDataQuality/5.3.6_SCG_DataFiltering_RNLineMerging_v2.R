#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)
IT = args[1]
print(IT)


library(ggplot2)
library(vroom)
library(tidyverse)
library(dplyr)
library(data.table)
library(circlize)
library(readxl)
library(dichromat)
library(RColorBrewer)
library(ggpubr)
library(stringr)
library(ComplexHeatmap)

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2")

#This loads and preps the data 
GROUPS <- vroom::vroom("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChrom_Model/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")
#GROUPS <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

MERGE_SCG_RPM <- vroom::vroom(paste("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/SingleCopyCoreGeneChunks/Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2_temp2_", IT, ".RPM.bed", sep = ""))

#This loads the data set with information on which lines passed missing data filters
PASS <- vroom::vroom("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/SamplesPASSEDMissingDataFILTER.csv")

#7 columns were classed as logicals because they have all NAs, these are the ones that had 0 coverage in the normalization process

#This drops any column that had an average coverage of 0, resulting in all values being NA
print("##### Starting Data Table Step")
DT <- as.data.table(MERGE_SCG_RPM)
DT <- DT[,which(unlist(lapply(DT, function(x)!all(is.na(x))))),with=F]
MERGE_SCG_RPM <- as.data.frame(DT)
print("##### Finished Data Table Step")

#This loads in the edited final BINS file
BINS <- read_excel("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/Bins_NoOverlap_SCG_final.xlsx")

print("##### Starting to add Bins")
#This adds the true bin to the SCG 
MERGE_SCG_RPM$bin <- NA
i=1
x=1
for (i in 1:nrow(MERGE_SCG_RPM)) {
  print(i)
  for (x in 1:nrow(BINS)) {
    BIN_NUM <- BINS[x,4][[1]]
    CHR <- BINS[x,1][[1]]
    START <- BINS[x,2][[1]]
    END <- BINS[x,3][[1]]
    if(MERGE_SCG_RPM[i,2]== CHR) {
      if (MERGE_SCG_RPM[i,3] > START & MERGE_SCG_RPM[i,4] < END) {
        MERGE_SCG_RPM[i,ncol(MERGE_SCG_RPM)] <- BIN_NUM
      }
    }
  }
}

print("##### Finsihed adding bins")
test <- MERGE_SCG_RPM[,c(1:6, ncol(MERGE_SCG_RPM))]

print("##### Filtering Low Coverage Samples")
#This selects only columns that passed the missing data filter
KEEP <- c(1:6, which(colnames(MERGE_SCG_RPM) %in% PASS$Name))
MERGE_SCG_RPM_FILT_1 <- MERGE_SCG_RPM[,KEEP]

print("##### Filtering Low MAPQ Tags")
#This removes any tag with a MAPQ less than 20
MERGE_SCG_RPM_FILT_2 <- subset(MERGE_SCG_RPM_FILT_1, MAPQ >= 20)

print("##### Merging the RN Samples")
#This generates a data frame with of all the line names found in the data set
NAME <- as.data.frame(colnames(MERGE_SCG_RPM_FILT_2[-c(1:7)]))
colnames(NAME) <- "Name"

#This loop alters the names of the Romero-Navarro data so that it can be merged with the groups file
NAME$Extra_Name <- NA
NAME$OG_Name<- NAME$Name
NAME$Simp_Name <- NA
i=1
for(i in 1:nrow(NAME)){
  SAMP <- NAME[i,1]
  if(endsWith(SAMP, "RN") == TRUE){
    TEMP <- str_split(SAMP, "[.]", n=2)
    NAME[i,4] <- TEMP[[1]][1]
    NAME[i,2] <- TEMP[[1]][2]
    NAME[i,3] <- SAMP
    NAME[i,1] <- paste((TEMP[[1]][1]),"RN", sep=".")
  }
}  

#This section isolates the Swarts Data and merges it with the groups.
GROUPS_SW <- subset(GROUPS, Data_Source == "Swarts_etal_2017")
NAME_SW <-merge(NAME, GROUPS_SW, by = "Name")
#this drops 27 samples which were part of the Swarts et al 2017 project, but not GBS sequenced

#This section isolates the Romay Data and merges it with the groups.
GROUPS_RY <- subset(GROUPS, Data_Source == "Romay_etal_2013")
NAME_RY <- merge(NAME, GROUPS_RY, by = "Name")

#This section isolates the Romay Data and merges it with the groups.
GROUPS_DC <- subset(GROUPS, Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2")
NAME_DC <- merge(NAME, GROUPS_DC, by = "Name")

#This section isolates the Romero Navarro Data and merges it with the groups.
GROUPS_RN <- subset(GROUPS, Data_Source == "Romero-Navarro_etal_2017")
NAME_RN <- merge((NAME[(grep("RN", NAME$Name)),]), GROUPS_RN, by.x = "Name", by.y="Name")

#This brings together all of the modified groups files with the altered names
DF <- rbind(NAME_SW, NAME_RY, NAME_RN, NAME_DC)

#This writes out the file for later use 
#fwrite(DF, file="Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table", row.names = FALSE)

##################This section isolates all of the data that does not need to be altered
#This identifies all the index numbers for non Romero-Navarro Lines

DF_OTH <- subset(DF, Data_Source != "Romero-Navarro_etal_2017")
col.num <- which(colnames(MERGE_SCG_RPM_FILT_2) %in% DF_OTH$OG_Name)
#This appends the meta data columns to this list
add <- c(1:6,ncol(MERGE_SCG_RPM_FILT_2))
col.num <- append(add, col.num)
#This selects all of these columns
MERGE_SCG_RPM_FILT_2_OTH <- MERGE_SCG_RPM_FILT_2[,c(col.num)]

###################This section takes the mean of all the rows for all the Romero-Navarro Data originating from the same individual

#This creates a list of all the unique names from the Romero-Navarro dataset.
SAMPS <- unique(NAME_RN$Name)

#This loop takes the mean of all normalized Romero-Navarro technical replicates
for(i in 1:length(SAMPS)) {
  x <- SAMPS[i]
  INDEX <- which(NAME_RN$Name %in% x)
  OGNAME <- NAME_RN[INDEX, "OG_Name"]
  INDEX2 <- which(colnames(MERGE_SCG_RPM_FILT_2) %in% OGNAME)
  TEMP1 <- as.data.frame(MERGE_SCG_RPM_FILT_2[,INDEX2], drop = FALSE)
  MEANS <- as.data.frame(rowMeans(TEMP1))
  CLEAN_NAME <- NAME_RN[INDEX, "Name"][1]
  colnames(MEANS) <- CLEAN_NAME
  MERGE_SCG_RPM_FILT_2_OTH <<- cbind(MERGE_SCG_RPM_FILT_2_OTH, MEANS)
}

MERGE_SCG_RPM_FILT_2_FIX <- MERGE_SCG_RPM_FILT_2_OTH

print("##### Writing the File")
#This writes out the final data frame
fwrite(MERGE_SCG_RPM_FILT_2_FIX, file=paste("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/SingleCopyCoreGeneChunks/BWAaln_All_v_Ab10HIFIBChrom.SCG.RPM.RNMean.", IT, ".table", sep = ""), row.names = FALSE)
