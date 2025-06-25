library("tidyverse")
library("readxl")

###### This Loads In Key Files
DC_1 <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/DaweLabControls1_Key.csv")

DC_2 <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/DaweLabControls2_Key.csv")

SW <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Swarts_Key.csv")

RN_key <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Romero-Navarro_LengthFiltered_fakebarcodes_key.csv")

Romay_list <- list.files(path="/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Romay_etal_2013_GBS/Keys")

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Romay_etal_2013_GBS/Keys")

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

##############This section makes the key for the Romay data

#Reads the whole file list into a single dataframe
Romay_Key_temp1 <- lapply(Romay_list,read_xlsx)
Romay_Key_temp2 <- dplyr::bind_rows(Romay_Key_temp1)
#Several rows are identical, I am not sure why this is. This removes any completely duplicate row
Romay_Key_temp3 <- unique(Romay_Key_temp2)

#This makes all the blank formating the same
Romay_Key_temp3$DNASample <- gsub("blank", "BLANK", Romay_Key_temp3$DNASample)
Romay_Key_temp3$DNASample <- gsub("Blank", "BLANK", Romay_Key_temp3$DNASample)

#This determiens the number of unique samepls
length(unique(Romay_Key_temp3$DNASample))

#This creates unique names for all Romay samples
SAMPS <- unique(Romay_Key_temp3$DNASample)
#This selects only the column names
Romay_Key <- Romay_Key_temp3[1,]
Romay_Key <- Romay_Key[-c(1),]

#This adds a number after duplicate samples to generate fully unique names
i=1
for(i in 1:length(SAMPS)) {
  x <- SAMPS[i]
  index <- which(Romay_Key_temp3$DNASample %in% x)
  TEMP1 <- as.data.frame(Romay_Key_temp3[index,], drop = FALSE)
  rownames(TEMP1) <- 1:nrow(TEMP1) 
  NEW <- paste(TEMP1$DNASample, rownames(TEMP1), sep=".")
  TEMP1$DNASample <- NEW
  Romay_Key <<- rbind(Romay_Key, TEMP1)
}

Romay_Key$LibraryPrepID <- Romay_Key$Flowcell
Romay_Key$FullSampleName <- Romay_Key$DNASample

#########################This section adds the data source initials to the end of each full sample name
DC_1$FullSampleName <- paste(DC_1$FullSampleName, "DC1", sep=".")
DC_2$FullSampleName <- paste(DC_2$FullSampleName, "DC2", sep=".")
SW$FullSampleName <- paste(SW$FullSampleName, "SW", sep=".")
RN_FinKey$FullSampleName <- paste(RN_FinKey$FullSampleName, "RN", sep=".")
Romay_Key$FullSampleName <- paste(Romay_Key$FullSampleName, "RY", sep=".")

##########################This section Combines the key files
All_Key <- rbind(DC_1, DC_2, SW, RN_FinKey, Romay_Key)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data")

write.table(x=All_Key, file="SwartsAllControlsLengthFiltRomeroNavarroRomay_Key.txt", row.names = FALSE, sep = "\t", quote = FALSE)


