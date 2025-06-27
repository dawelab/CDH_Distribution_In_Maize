#Load the finalized file
DATA <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumComplete.csv")

#Extract Ab10 samples only
DATA_Ab10 <- subset(DATA, KMeans_Ab10 == "Ab10")
DATA_Ab10 <- DATA_Ab10[c("KMeans_Ab10", "Longitude", "Latitude")]
names(DATA_Ab10) <- c("SGE", "Longitude", "Latitude")

#Extract K10L2 samples only
DATA_K10L2 <- subset(DATA, KMeans_K10L2 == "K10L2")
DATA_K10L2 <- DATA_K10L2[c("KMeans_K10L2", "Longitude", "Latitude")]
names(DATA_K10L2) <- c("SGE", "Longitude", "Latitude")

#Extract B Chrom samples
DATA_B <- subset(DATA, KMeans_BChrom == "Yes")
DATA_B$KMeans_BChrom <- "BChrom" 
DATA_B <- DATA_B[c("KMeans_BChrom", "Longitude", "Latitude")]
names(DATA_B) <- c("SGE", "Longitude", "Latitude")

#Extract GPS for all samples
DATA_ALL <- DATA[c("Longitude", "Latitude")]
DATA_ALL$SGE <- "Maize"

#Combine all the files
DATA_ALL <- rbind(DATA_Ab10, DATA_K10L2, DATA_B, DATA_ALL)

#This drops any sample without GPS coordinates
DATA_ALL_FILT <- subset(DATA_ALL, is.na(DATA_ALL$Longitude) == FALSE & is.na(DATA_ALL$Latitude) == FALSE )

#Write out the formatted file 
write.csv(DATA_ALL_FILT, "SGE_GPS_MaxEnt_Format.csv", row.names = FALSE)
