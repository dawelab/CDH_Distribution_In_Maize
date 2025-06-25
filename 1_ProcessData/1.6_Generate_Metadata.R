library(dplyr)

#Load the relevant files
MAP <- read_excel("Mapping file - Sample ID to Germplasm ID download2.xlsx")
GERM <- read_excel("Germinate_Group5_GbS_Accessions_MBEdit.xlsx")
file_num <- read.table("file_num.txt", quote="\"", comment.char="")
names(file_num) <- c("Permissions", "Nothing", "User", "PI", "Size", "Month", "Day", "Time", "Sample ID")
Swarts_Groups_Env <- as.data.frame(read.csv("SwartsControls_Groups_Env.csv"))


#There are multiple instances of general identifiers in the mapping file. 
length(unique(MAP$general_identifier))
length(unique(GERM$general_identifier))


length(MAP$general_identifier) - length(unique(MAP$general_identifier))
#690 samples are replicates of the same general identifier, I believe that these are because they were used to make different testcrosses for their FOAM population

#This merges the two data frames
TEMP <- merge(MAP, GERM, by = "general_identifier")

#This extracts only lines that also have GBS data
ALL <- left_join(file_num, TEMP, by = "Sample ID")
length(unique(ALL$general_identifier))

GPS_only <- subset(ALL, is.na(ALL$locations_latitude)==FALSE)
length(ALL$general_identifier)-length(GPS_only$general_identifier)

length(unique(GPS_only$general_identifier))
length(ALL$general_identifier)-length(GPS_only$general_identifier)

#There are 4021 unqiue samples with Germinate data and GPS data. Of these 2898 have GPS coordinates associated. 1123 do not have GPS coordinates

#There are 3523 lines that have both GBS data and GPS coordinates, of these 2898 have unique general identifiers. There are 1190 samples that have GBS data but do not have GPS coordinates. 

#I now need to put these into the same format at the other data
Name <- ALL$`Sample ID`
Data_Source <- "Romero-Navarro"
Accession <- ALL$bank_number
Group <- "NA"
Maize_Type <- "Landrace"
Ab10_Status <- "Unknown"
Ab10_Status_Spec <- "Unknown"
Latitude <- ALL$locations_latitude
Longitude <- ALL$locations_longitude
Altitude <- ALL$locations_elevation
DTA_BLUP <- "NA"
DTS_BLUP <- "NA"
PH_BLUP <- "NA"
PresTiller_BLUP <- "NA"

BOUND <- cbind(Name, Data_Source, Accession, Group, Maize_Type, Ab10_Status, Ab10_Status_Spec, Latitude, Longitude, Altitude, DTA_BLUP, DTS_BLUP, PH_BLUP, PresTiller_BLUP)
RN_DF <- as.data.frame(BOUND)

Swarts_RomeroNavarro_Groups_Env <- rbind(Swarts_Groups_Env, RN_DF)

write.csv(Swarts_RomeroNavarro_Groups_Env, file="Swarts_Controls_RomeroNavarro_Groups_Env.csv")
