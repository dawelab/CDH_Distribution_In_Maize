library(stringr)

SGE <- read.csv("/Volumes/Transcend/6.5_Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumComplete_ForSGE.csv")
KEY <- read.delim("/Volumes/Transcend/SwartsAllControlsLengthFiltRomeroNavarroRomay_Key_ForSGE.txt")
AB10 <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/All_Ab10_Positive.txt", sep="")

SGE <- subset(SGE, SGE$KMeans_Ab10 == "Ab10" | SGE$KMeans_K10L2 == "K10L2" | SGE$KMeans_BChrom == "Yes")

SGE_KEY <- KEY[1,]
SGE_KEY <- SGE_KEY[-c(1),]

for(i in 1:nrow(SGE)) {
  x <- SGE[i,1]
  NAME <- stringr::str_split(x, "[.]")[[1]][1]
  SOURCE <- stringr::str_split(x, "[.]")[[1]][2]
  SUB_TEMP <- subset(KEY, grepl(NAME, KEY$FullSampleName))
  SUB <- subset(SUB_TEMP, grepl(SOURCE, SUB_TEMP$FullSampleName))
  SUB$FullSampleName <- paste(NAME, SOURCE, sep = ".")
  SGE_KEY <<- rbind(SGE_KEY, SUB)
}

length(unique(SGE_KEY$FullSampleName))

write.table(SGE_KEY, file="/Volumes/Transcend/SwartsAllControlsLengthFiltRomeroNavarroRomay_Key_SGEOnly.txt", quote = FALSE, row.names = FALSE)

