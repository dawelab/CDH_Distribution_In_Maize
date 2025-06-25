#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)
DIRNAME = args[1]
print(DIRNAME)
FILNAME = args[2]
print(FILNAME)

FILNAME<-"Tassel_TagTaxaDist_AllData_v2_v_B73-Ab10_BChrom.Sub1Perc.15.txt"

library(vroom)
library(data.table)

setwd(DIRNAME)

fil <- vroom::vroom(FILNAME)

#This drops the first column containing tag information 
fil <- fil[,-c(1)]

#This saves the column names as a vector
Names <- colnames(fil)
#This calculates sums across
Sums <- apply(fil, 2, sum)

SUM_DATA <- as.data.frame(cbind(Names, Sums))

SHORTNAME <- strsplit(FILNAME, "[.]")[[1]][1]
IT <- strsplit(FILNAME, "[.]")[[1]][3]

fwrite(SUM_DATA, file=paste(SHORTNAME, ".Sub1Perc.Sum.", IT, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
