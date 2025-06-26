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

#This loads and preps the data 
#This file is from 1.7
GROUPS <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

SCG_GFF <- vroom::vroom("Ab10_HiFi_v2_corrected.liftoff.CoreSingleCopy.gff3")
colnames(SCG_GFF) <- c("chr", "Liftoff", "Feature", "start", "end", "Dot1", "Strand", "Dot2", "ID")

SCG_GFF$length <- SCG_GFF$end - SCG_GFF$start
sum(SCG_GFF$length)

SCG_GFF$pseudo_start <- "NA"
SCG_GFF$pseudo_end <- "NA"

write.csv(SCG_GFF, file = "Ab10_HiFi_v2_corrected.liftoff.CoreSingleCopy.temp1.gff3.csv", row.names = FALSE, quote = FALSE)

#I went into excel and created the fields pseudo start and end 
SCG_GFF <- vroom::vroom("Ab10_HiFi_v2_corrected.liftoff.CoreSingleCopy.temp2.gff3.csv")

#This loads in the BINS file
BINS <- read_excel("Bins_NoOverlap_SCG.xlsx")

#this adds a bin value to each line in the SCG_GFF genes 
SCG_GFF$bin <- NA
i=1
x=1
for (i in 1:nrow(SCG_GFF)) {
  print(i)
  for (x in 1:nrow(BINS)) {
    BIN_NUM <- BINS[x,4][[1]]
    START <- BINS[x,2][[1]]
    END <- BINS[x,3][[1]]
    if (SCG_GFF[i,11] > START & SCG_GFF[i,12] < END) {
      SCG_GFF[i,ncol(SCG_GFF)] <- BIN_NUM
    }
  }
}

#This writes out the file
write.csv(SCG_GFF, file ="Ab10_HiFi_v2_corrected.liftoff.CoreSingleCopy.bin.gff3.csv", row.names = FALSE, quote = FALSE)

#I went in and generated a bins file from the one above which ensures that there is 1MB of single copy core gene genomic sequence in each bin regardless of their physical distance

#This loads in the edited final BINS file, this file is available in this repo under 5.4.5
BINS <- read_excel("Bins_NoOverlap_SCG_final.xlsx")
