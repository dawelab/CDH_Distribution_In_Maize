library("tidyverse")
library("readxl")

###### This Loads In Key Files

RN_key <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Romero-Navarro_LengthFiltered_fakebarcodes_key.csv")

######################This section makes the key for the RomeroNavarro data 
Flowcell <- RN_key$fake_flowcells
Lane <- RN_key$fake_lanes
Barcode <- RN_key$fake_barcodes
DNAsample_temp1 <- RN_key$fastqs
DNASample <- c()
FullSampleName <- c()
LibraryPrepID <- 1:length(DNAsample_temp1)

for (i in 1:length(DNAsample_temp1)) {  
SAMP <- DNAsample_temp1[[i]]
CUT1 <- str_split(SAMP, "/")[[1]][7]
CUT2 <- str_split(CUT1, "_")[[1]][2]
CUT3 <- str_split(CUT2, ":")[[1]][1]
DNASample[i] <- CUT2
FullSampleName[i] <- CUT3
}

RN_NewKey <- cbind(Flowcell, Lane, Barcode, DNASample, LibraryPrepID, FullSampleName)
RN_NewKey <- as.data.frame(RN_NewKey)

#This creates unique names for all Romay samples
SAMPS <- unique(RN_NewKey$FullSampleName)
#This selects only the column names
RN_FinKey <- RN_NewKey[1,]
RN_FinKey <- RN_FinKey[-c(1),]

#This adds a number after duplicate samples to generate fully unique names
i=1
for(i in 1:length(SAMPS)) {
  x <- SAMPS[i]
  index <- which(RN_NewKey$FullSampleName %in% x)
  TEMP1 <- as.data.frame(RN_NewKey[index,], drop = FALSE)
  rownames(TEMP1) <- 1:nrow(TEMP1) 
  NEW <- paste(TEMP1$FullSampleName, rownames(TEMP1), sep=".")
  TEMP1$FullSampleName <- NEW
  RN_FinKey <<- rbind(RN_FinKey, TEMP1)
}

#This writes out a file to fix the Romero Navarro names. for elsewhere. It is only needed because of an error in the key
RN_FinKey_Fix <- RN_FinKey
RN_FinKey_Fix$FullSampleName_Fix <- paste(RN_FinKey$FullSampleName, "RN", sep=".")
setwd("/Volumes/Transcend/")
write.csv(RN_FinKey_Fix, file="Romero_Navarro_Key_Fix.csv", row.names = FALSE, quote = FALSE)
