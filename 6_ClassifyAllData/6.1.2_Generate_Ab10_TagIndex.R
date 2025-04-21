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
setwd("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2")

#This loads the original groups file
GROUPS <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This loads in the filtered, Romero Navarro merged data set. 
MERGE_Ab10Hap_RPM_FILT_2_FIX <- vroom::vroom("Ab10_Model/BWAaln_All_v_Ab10HIFIBChrom.Ab10.RPM.RNMean.table")

#This loads in the groups file with modified names
DF <- vroom::vroom("Ab10_Model/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table")

#This loads in the BINS file
BINS <- read_excel("Bins_NoOverlap.table.xlsx")

#This created output directories
#dir.create("Ab10_Model/Experimental_Ab10N10Model/Images")
#dir.create("Ab10_Model/Experimental_Ab10N10Model/SubSamples")

#This function goes over each row and divides each value by the max in that row 
MinMax = function(xx) { sweep(xx, 1, apply(xx, 1, max), '/') }

######### Tag Sums
#This sums the number of tag numbers across bins, but excludes the bin column from the summing

A_SUM = as.data.frame(apply(MERGE_Ab10Hap_RPM_FILT_2_FIX[,-c(1:7)], 2, function(xx) { by(xx, MERGE_Ab10Hap_RPM_FILT_2_FIX$bin, sum)  }))

#This writes out the sumed file
fwrite(A_SUM, file ="Ab10_Model/Experimental_Ab10N10Model/Ab10_TagSums_AllSamples.csv")

######## Tag Density
#This sums the number of tag numbers across bins, but excludes the bin column from the summing

A_DENS = as.data.frame(apply(MERGE_Ab10Hap_RPM_FILT_2_FIX[,-c(1:7)], 2, function(xx) { by(xx, MERGE_Ab10Hap_RPM_FILT_2_FIX$bin, function(yy) { sum(yy > 0) }) }))

#This writes out the tag density file
fwrite(A_DENS, file ="Ab10_Model/Experimental_Ab10N10Model/Ab10_TagDens_AllSamples.csv")

#This takes the square root of each value
B_DENS <- sqrt(A_DENS)

#This combines the density and sum data
B_COMB <- A_SUM+B_DENS

fwrite(B_COMB, file ="Ab10_Model/Experimental_Ab10N10Model/Ab10_TagSumPlusTagDensAllSamples.csv")
