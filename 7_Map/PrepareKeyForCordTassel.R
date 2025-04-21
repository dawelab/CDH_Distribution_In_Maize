library(stringr)
library(dplyr)

SGE <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/6.5_Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumComplete.csv")
KEY <- read.delim("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/SwartsAllControlsLengthFiltRomeroNavarroRomay_Key.txt")

#For some reason the key has duplicate records. I think it was corrupted in download and processing. I am removing duplicates
KEY <- distinct(KEY)

#This subsets to only the lines with GPS coordinates
CORD <- subset(SGE, is.na(SGE$Latitude) == FALSE | is.na(SGE$Longitude) == FALSE)

#This preps the new dataframe
CORD_KEY <- CORD[1,]
CORD_KEY <- CORD_KEY[-c(1),]

#This loop edits the names so that all sequences from the same biological sample have the same full name
for(i in 1:nrow(CORD)) {
  x <- CORD[i,1]
  NAME <- stringr::str_split(x, "[.]")[[1]][1]
  SOURCE <- stringr::str_split(x, "[.]")[[1]][2]
  SUB_TEMP <- subset(KEY, grepl(paste(NAME, ".", sep=""), KEY$FullSampleName, fixed = TRUE))
  SUB <- subset(SUB_TEMP, grepl(SOURCE, SUB_TEMP$FullSampleName))
  SUB$FullSampleName <- paste(NAME, SOURCE, sep = ".")
  CORD_KEY <<- rbind(CORD_KEY, SUB)
}

length(unique(CORD_KEY$FullSampleName))

write.table(CORD_KEY, file="/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/SwartsAllControlsLengthFiltRomeroNavarroRomay_Key_CordSNP.txt", quote = FALSE, row.names = FALSE)
