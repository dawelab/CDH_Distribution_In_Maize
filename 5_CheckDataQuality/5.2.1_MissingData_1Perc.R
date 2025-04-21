#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)
DIRNAME = args[1]
print(DIRNAME)
FILNAME = args[2]
print(FILNAME)

library(vroom)
library(data.table)

setwd(DIRNAME)

fil <- vroom::vroom(FILNAME)

#This drops the first column containing tag information 
fil <- fil[,-c(1)]

#This saves the column names as a vector
Names <- colnames(fil)
#This calculates missing data
Miss <- colSums(fil == 0)
Total <- nrow(fil)
Perc_Miss <- Miss/Total

MISS_DATA <- as.data.frame(cbind(Names, Miss, Total, Perc_Miss))

SHORTNAME <- strsplit(FILNAME, "[.]")[[1]][1]
IT <- strsplit(FILNAME, "[.]")[[1]][3]

fwrite(MISS_DATA, file=paste(SHORTNAME, ".Sub1Perc.Missing.", IT, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
