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

#This sets the working directory where all files being loaded in are located
setwd("")

#This loads and preps the data 
#This file is from 1.7
GROUPS <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This file is from 3.4
MERGE_Ab10Hap_RPM <- vroom::vroom("Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10_BChrom.Ab10.RPM.txt")
#7 columns were classed as logicals because they have all NAs, these are the ones that had 0 coverage

#This loads the data set with information on which lines passed missing data filters
#This file is from 5.3
PASS <- vroom::vroom("SamplesPASSEDMissingDataFILTER.csv")

#This drops any column that had an average coverage of 0, resulting in all values being NA
DT <- as.data.table(MERGE_Ab10Hap_RPM)
DT <- DT[,which(unlist(lapply(DT, function(x)!all(is.na(x))))),with=F]
MERGE_Ab10Hap_RPM <- as.data.frame(DT)

#This loads in the BINS file
#This file is available in 4.1
BINS <- read_excel("Bins_NoOverlap.table.xlsx")

#This selects only columns that passed the missing data filter
KEEP <- c(1:6, which(colnames(MERGE_Ab10Hap_RPM) %in% PASS$Name))
MERGE_Ab10Hap_RPM_FILT_1 <- MERGE_Ab10Hap_RPM[,KEEP]

#This removes any tag with a MAPQ less than 20
MERGE_Ab10Hap_RPM_FILT_2 <- subset(MERGE_Ab10Hap_RPM_FILT_1, MAPQ >= 20)

#this adds a bin value to each line in the MERGE_Ab10Hap_RPM file 
MERGE_Ab10Hap_RPM_FILT_2$bin <- NA
i=1
x=1
for (i in 1:nrow(MERGE_Ab10Hap_RPM_FILT_2)) {
  print(i)
  for (x in 1:nrow(BINS)) {
    BIN_NUM <- BINS[x,4][[1]]
    START <- BINS[x,2][[1]]
    END <- BINS[x,3][[1]]
    if (MERGE_Ab10Hap_RPM_FILT_2[i,3] > START & MERGE_Ab10Hap_RPM_FILT_2[i,3] < END) {
      MERGE_Ab10Hap_RPM_FILT_2[i,ncol(MERGE_Ab10Hap_RPM_FILT_2)] <- BIN_NUM
    }
  }
}

#This generates a data frame with of all the line names found in the data sest
NAME <- as.data.frame(colnames(MERGE_Ab10Hap_RPM_FILT_2[-c(1:7)]))
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
fwrite(DF, file="Ab10_Model/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table", row.names = FALSE)

##################This section isolates all of the data that does not need to be altered
#This identifies all the index numbers for non Romero-Navarro Lines

DF_OTH <- subset(DF, Data_Source != "Romero-Navarro_etal_2017")
col.num <- which(colnames(MERGE_Ab10Hap_RPM_FILT_2) %in% DF_OTH$OG_Name)
#This appends the meta data columns to this list
add <- c(1:6,ncol(MERGE_Ab10Hap_RPM_FILT_2))
col.num <- append(add, col.num)
#This selects all of these columns
MERGE_Ab10Hap_RPM_FILT_2_OTH <- MERGE_Ab10Hap_RPM_FILT_2[,c(col.num)]


###################This section takes the mean of all the rows for all the Romero-Navarro Data originating from the same individual

#This creates a list of all the unique names from the Romero-Navarro dataset.
SAMPS <- unique(NAME_RN$Name)

#This loop takes the mean of all normalized Romero-Navarro technical replicates
for(i in 1:length(SAMPS)) {
  x <- SAMPS[i]
  INDEX <- which(NAME_RN$Name %in% x)
  OGNAME <- NAME_RN[INDEX, "OG_Name"]
  INDEX2 <- which(colnames(MERGE_Ab10Hap_RPM_FILT_2) %in% OGNAME)
  TEMP1 <- as.data.frame(MERGE_Ab10Hap_RPM_FILT_2[,INDEX2], drop = FALSE)
  MEANS <- as.data.frame(rowMeans(TEMP1))
  CLEAN_NAME <- NAME_RN[INDEX, "Name"][1]
  colnames(MEANS) <- CLEAN_NAME
  MERGE_Ab10Hap_RPM_FILT_2_OTH <<- cbind(MERGE_Ab10Hap_RPM_FILT_2_OTH, MEANS)
}

MERGE_Ab10Hap_RPM_FILT_2_FIX <- MERGE_Ab10Hap_RPM_FILT_2_OTH

#This writes out the final data frame
fwrite(MERGE_Ab10Hap_RPM_FILT_2_FIX, file="Ab10_Model/BWAaln_All_v_B73-Ab10_BChrom.Ab10.RPM.RNMean.table", row.names = FALSE)
