library( StAMPP)
library(stringr)

#I manually removed the header from the file output by Filter_Higgins_SNPs.R

VCF <- read.table("/Volumes/Transcend/AllData_v5.Ab10Hap.HigControls.filt.nohead.vcf", header=TRUE)

#This function converts the filtered vcf to AB format with -9 as the missing symbol for compatibillity with StAMPP
convert <- function(x) {
  y <- str_split((x), ":")[[1]][1]
  y <- sub("0/1", "AB", y)
  y <- sub("1/1", "BB", y)
  y <- sub("0/0", "AA", y)
  y <- sub("1/0", "AB", y)
  y <- sub("./.", -9, y)
  x <<-y
}

VCF[,-c(1:9)] <- apply(VCF[,-c(1:9)], c(1, 2), convert)

ID <- VCF$ID

#This reformats the data for compatibility with StAMPP
VCF <- VCF[,-c(1:9)]
#This transposes the data frame
VCF <- (t(VCF))
#This converts the VCF to a data frame
VCF <- as.data.frame(VCF)
#This assigns columns names
colnames(VCF) <- c(1:ncol(VCF))
#This assigns row names to a new field
VCF$Sample <- rownames(VCF)
#This assigns a population value
VCF$Pop <- c(rep("Ab10-1", 10), rep("Ab10-2", 15), rep("Ab10-3", 20))
#This assigns a ploidy value
VCF$Ploidy <- 2
#This assigns the format
VCF$Format <- "	BiA"
#This reorders the data frame
VCF <- VCF[,c("Sample", "Pop", "Ploidy", "Format", 1:591)]
#This renames the SNPs to reflect their position
colnames(VCF) <- c("Sample", "Pop", "Ploidy", "Format", ID)

#This writes out the file
write.table(VCF, file="/Volumes/Transcend/AllData_v5.Ab10Hap.HigControls.filt.AB.vcf", quote=FALSE, row.names = FALSE)


