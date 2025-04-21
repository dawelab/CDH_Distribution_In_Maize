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

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/")

#This loads the original groups file
GROUPS <- vroom::vroom("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChrom_Model/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This reads in all the filtered files together  
IT=1
MERGE_SCG_RPM_FILT_2_FIX <- vroom::vroom(paste("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/SingleCopyCoreGeneChunks/BWAaln_All_v_Ab10HIFIBChrom.SCG.RPM.RNMean.", IT, ".table", sep = ""))

for(i in 1:19) {
  i=i+1
  print(i)
  TEMP <- vroom::vroom(paste("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/SingleCopyCoreGeneChunks/BWAaln_All_v_Ab10HIFIBChrom.SCG.RPM.RNMean.", i, ".table", sep = ""))
  MERGE_SCG_RPM_FILT_2_FIX <<- rbind(MERGE_SCG_RPM_FILT_2_FIX, TEMP)
}

#This writes out the full file
fwrite(MERGE_SCG_RPM_FILT_2_FIX, file="BWAaln_All_v_Ab10HIFIBChrom.SCG.RPM.RNMean.All.table")

#This loads in the groups file with modified names
DF <- vroom::vroom("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChrom_Model/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table")

#This loads in the BINS file
BINS <- read_excel("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChrom_Model/Bins_NoOverlap_BChrom.table.xlsx")

#This function goes over each row and divides each value by the max in that row 
MinMax = function(xx) { sweep(xx, 1, apply(xx, 1, max), '/') }

######### Tag Sums
#This sums the number of tag numbers across bins, but excludes the bin column from the summing

A_SUM = as.data.frame(apply(MERGE_SCG_RPM_FILT_2_FIX[,-c(1:7)], 2, function(xx) { by(xx, MERGE_SCG_RPM_FILT_2_FIX$bin, sum)  }))

#This writes out the sumed file
fwrite(A_SUM, file ="SCG_TagSums_AllSamples.csv")

######## Tag Density
#This sums the number of tag numbers across bins, but excludes the bin column from the summing

A_DENS = as.data.frame(apply(MERGE_SCG_RPM_FILT_2_FIX[,-c(1:7)], 2, function(xx) { by(xx, MERGE_SCG_RPM_FILT_2_FIX$bin, function(yy) { sum(yy > 0) }) }))


#This writes out the tag density file
fwrite(A_DENS, file ="SCG_TagDens_AllSamples.csv")

#This takes the square root of each value
B_DENS <- sqrt(A_DENS)

#This combines the density and sum data
B_COMB <- A_SUM+B_DENS

fwrite(B_COMB, file ="SCG_TagSumPlusTagDensAllSamples.csv")
