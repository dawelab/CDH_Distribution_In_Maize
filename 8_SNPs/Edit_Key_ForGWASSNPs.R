library(stringr)

GROUPS <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/6.5_Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumComplete.csv")
KEY <- read.delim("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/SwartsAllControlsLengthFiltRomeroNavarroRomay_Key.txt")

#This drops duplicate rows
KEY <- KEY[!duplicated(KEY), ]

#This preps an empty key
GWAS_KEY <- KEY[1,]
GWAS_KEY <- GWAS_KEY[-c(1),]

#This edits the Romay samples in the GROUPS file
RY <- GROUPS[grep("RY", GROUPS$Name),]
GROUPS <- GROUPS[-grep("RY", GROUPS$Name),]

for(i in 1:nrow(RY)) {
a <- RY[i,1]
NAME <- stringr::str_split(a, "[.]")[[1]][1]
NUM <- stringr::str_split(a, "[.]")[[1]][2]
SOURCE <- stringr::str_split(a, "[.]")[[1]][3]
RY[i,1] <- paste(NAME, "_", NUM, ".", SOURCE, sep = "")
}

GROUPS <- rbind(GROUPS, RY)

#This edits the Romay samples in the Key file
RY <- KEY[grep("RY", KEY$FullSampleName),]
KEY <- KEY[-grep("RY", KEY$FullSampleName),]

for(i in 1:nrow(RY)) {
  a <- RY[i,6]
  NAME <- stringr::str_split(a, "[.]")[[1]][1]
  NUM <- stringr::str_split(a, "[.]")[[1]][2]
  SOURCE <- stringr::str_split(a, "[.]")[[1]][3]
  RY[i,6] <- paste(NAME, "_", NUM, ".", SOURCE, sep = "")
}

KEY <- rbind(KEY, RY)

#This edits the RIMMA samples in the GROUPS file
MA <- GROUPS[grep("RIMMA", GROUPS$Name),]
GROUPS <- GROUPS[-grep("RIMMA", GROUPS$Name),]

for(i in 1:nrow(MA)) {
  a <- MA[i,1]
  NAME <- stringr::str_split(a, "[.]")[[1]][1]
  NUM <- stringr::str_split(a, "[.]")[[1]][2]
  SOURCE <- stringr::str_split(a, "[.]")[[1]][3]
  MA[i,1] <- paste(NAME, "_", NUM, ".", SOURCE, sep = "")
}

GROUPS <- rbind(GROUPS, MA)

#This edits the RIMMA samples in the Key file
MA <- KEY[grep("RIMMA", KEY$FullSampleName),]
KEY <- KEY[-grep("RIMMA", KEY$FullSampleName),]

for(i in 1:nrow(MA)) {
  a <- MA[i,6]
  NAME <- stringr::str_split(a, "[.]")[[1]][1]
  NUM <- stringr::str_split(a, "[.]")[[1]][2]
  SOURCE <- stringr::str_split(a, "[.]")[[1]][3]
  MA[i,6] <- paste(NAME, "_", NUM, ".", SOURCE, sep = "")
}

KEY <- rbind(KEY, MA)

#This merges duplicate samples for RN
MA <- KEY[grep(".RN", KEY$FullSampleName),]
GWAS_KEY <- KEY[-grep(".RN", KEY$FullSampleName),]

for(i in 1:nrow(MA)) {
  x <- MA[i,6]
  NAME <- stringr::str_split(x, "[.]")[[1]][1]
  SOURCE <- stringr::str_split(x, "[.]")[[1]][3]
  SUB_TEMP <- subset(MA, grepl(paste0(NAME,"\\."), MA$FullSampleName))
  SUB <- subset(SUB_TEMP, grepl(SOURCE, SUB_TEMP$FullSampleName))
  SUB$FullSampleName <- paste(NAME, SOURCE, sep = ".")
  GWAS_KEY <<- rbind(GWAS_KEY, SUB)
}


#This drops duplicate rows
GWAS_KEY <- GWAS_KEY[!duplicated(GWAS_KEY), ]

write.table(GWAS_KEY, file="/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/8_SNPs/SwartsAllControlsLengthFiltRomeroNavarroRomay_Key_GWASSNPs.txt", quote = FALSE, row.names = FALSE, sep = "\t")

