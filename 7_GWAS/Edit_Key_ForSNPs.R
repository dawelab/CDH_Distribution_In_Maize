library(stringr)

SGE <- read.csv("/Volumes/Transcend/6.5_Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumComplete_ForSGE.csv")
KEY <- read.delim("/Volumes/Transcend/SwartsAllControlsLengthFiltRomeroNavarroRomay_Key_ForSGE.txt")
AB10 <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/All_Ab10_Positive.txt", sep="")

SGE <- subset(SGE, SGE$KMeans_Ab10 == "Ab10" | SGE$KMeans_K10L2 == "K10L2" | SGE$KMeans_BChrom == "Yes")

SGE_KEY <- KEY[1,]
SGE_KEY <- SGE_KEY[-c(1),]

#This edits the Romay samples in the SGE file
RY <- SGE[grep("RY", SGE$Name),]
SGE <- SGE[-grep("RY", SGE$Name),]

for(i in 1:nrow(RY)) {
a <- RY[i,1]
NAME <- stringr::str_split(a, "[.]")[[1]][1]
NUM <- stringr::str_split(a, "[.]")[[1]][2]
SOURCE <- stringr::str_split(a, "[.]")[[1]][3]
RY[i,1] <- paste(NAME, "_", NUM, ".", SOURCE, sep = "")
}

SGE <- rbind(SGE, RY)

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

#This edits the RIMMA samples in the SGE file
MA <- SGE[grep("RIMMA", SGE$Name),]
SGE <- SGE[-grep("RIMMA", SGE$Name),]

for(i in 1:nrow(MA)) {
  a <- MA[i,1]
  NAME <- stringr::str_split(a, "[.]")[[1]][1]
  NUM <- stringr::str_split(a, "[.]")[[1]][2]
  SOURCE <- stringr::str_split(a, "[.]")[[1]][3]
  MA[i,1] <- paste(NAME, "_", NUM, ".", SOURCE, sep = "")
}

SGE <- rbind(SGE, MA)

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

#This merges duplicate samples
for(i in 1:nrow(SGE)) {
  x <- SGE[i,1]
  NAME <- stringr::str_split(x, "[.]")[[1]][1]
  SOURCE <- stringr::str_split(x, "[.]")[[1]][2]
  SUB_TEMP <- subset(KEY, grepl(NAME, KEY$FullSampleName))
  SUB <- subset(SUB_TEMP, grepl(SOURCE, SUB_TEMP$FullSampleName))
  SUB$FullSampleName <- paste(NAME, SOURCE, sep = ".")
  SGE_KEY <<- rbind(SGE_KEY, SUB)
}

SGE_Ab10 <- subset(SGE, SGE$KMeans_Ab10 == "Ab10")
write.table(SGE_Ab10$Name, file="/Volumes/Transcend/All_Ab10_Pos.txt", quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)


SGE_K10L2 <- subset(SGE, SGE$KMeans_K10L2 == "K10L2")
write.table(SGE_K10L2$Name, file="/Volumes/Transcend/All_K10L2_Pos.txt", quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)

write.table(SGE_KEY, file="/Volumes/Transcend/SwartsAllControlsLengthFiltRomeroNavarroRomay_Key_SGEOnly.txt", quote = FALSE, row.names = FALSE, sep = "\t")

