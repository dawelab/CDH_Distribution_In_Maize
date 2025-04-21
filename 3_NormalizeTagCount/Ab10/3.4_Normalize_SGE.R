#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)
DIRNAME <-args[1]
print(DIRNAME)
NAME = args[2]
print(NAME)
SGE_N = args[3]
print(SGE_N)

library(vroom)
library(tidyverse)
library(data.table)

setwd(DIRNAME)

SGE <- vroom::vroom(paste(NAME,SGE_N,"txt", sep="."))

BEDNAME = gsub("Tassel_TagTaxaDist", "BWAaln", NAME)

BED <- vroom::vroom(paste(BEDNAME, "s.bed", sep="."))
colnames(BED) <- c("Chr", "Start", "End", "Tag", "MAPQ", "Strand")

BED$Tag <- gsub("tagSeq=", "", BED$Tag)

MERGE <- merge(BED, SGE, by="Tag")

###################This section of code selects only tags that appear in the control data set 

#This selects columns ending in the control specific suffix
col.num <- grep("DC", colnames(MERGE))
#This appends the first 6 columns to that list
col.num <- c(1:6, col.num)
#This selects only the columns that appear in the above list of control names
MERGE_CONTROLS <- MERGE[,sort(c(col.num))]
#This calculates the row sum for each tag
MERGE_CONTROLS$RowSum <- rowSums(MERGE_CONTROLS[7:ncol(MERGE_CONTROLS)])
#This removes any column where all of the tag values are 0 in the controls
MERGE_CONTROLS_NOZERO <- subset(MERGE_CONTROLS, RowSum != 0)
#This stores the tags present in the control data as a variable
TAGS_TO_KEEP <- MERGE_CONTROLS_NOZERO$Tag

#This goes back to the full dataset and stores tag information as row names
row.names(MERGE) <- MERGE$Tag
#This goes back to the original data and pulls the index for only tags that appear in the controls 
row.num <- which(rownames(MERGE) %in% TAGS_TO_KEEP)
#This selects only the columns that appear in the above list of tag index
MERGE_TAGSFilt <- MERGE[sort(c(row.num)),]

#This section loads in the sum data file containing the sums of 1% of a random subset of the data 

SUM_DATA <- vroom::vroom(paste(NAME, "Sub1Perc.Sum.txt", sep="."))

SUM_DATA$Mean = apply(SUM_DATA[,2:ncol(SUM_DATA)], 1, mean)

fwrite(SUM_DATA, file=paste(NAME, "Sub1Perc.Sum.Mean.txt", sep="."), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

###################This section normalizes tag count by the total tags. I tested a few ways to do this and this seemed to be the fastest

COLNAMES <- colnames(MERGE_TAGSFilt)
MERGE_TAGSFilt_Norm <- MERGE_TAGSFilt[1:6]

c<- 1
normalize <- function(x) {
  print(c)
  NAME <- x
  SUMROW <- which(SUM_DATA$Names == NAME)
  SUM <- as.numeric(SUM_DATA[SUMROW, ncol(SUM_DATA)])*100
  #This takes a vector and normalizes across rows 
  z <- which(COLNAMES == NAME)
  SUB <- as.data.frame(MERGE_TAGSFilt[,z:z])
  edit <- apply(SUB, 1, function(y) (y/SUM)*1000000)
  MERGE_TAGSFilt_Norm <<- cbind(MERGE_TAGSFilt_Norm, edit)
  c<<-c+1
}

COLNAMES_Drop <- COLNAMES[7:length(COLNAMES)]
lapply(COLNAMES_Drop, normalize)

#This assigns the colnames back to the data frame after transformation 
colnames(MERGE_TAGSFilt_Norm) <- COLNAMES

fwrite(MERGE_TAGSFilt_Norm, file=paste(NAME,SGE_N,"RPM", "txt", sep="."), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
