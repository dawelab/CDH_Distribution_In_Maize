library(ggplot2)

#This loads in the joined files with the Missing Data
MISS <- read.table("Tassel_TagTaxaDist_AllData_v7_v_K10L2.Sub1Perc.Missing.txt", header = TRUE)

#This loads the Groups File
#This file is from 1.7
GROUPS <- read.csv("Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This isolates the Names
MISS_v2 <- MISS[,c("Names")]

#This isolates and sums the Missing data count
MISS_v3 <- MISS[,c("Miss", "Miss.1", "Miss.2", "Miss.3", "Miss.4", "Miss.5", "Miss.6", "Miss.7", "Miss.8", "Miss.9", "Miss.10", "Miss.11", "Miss.12", "Miss.13", "Miss.14")]
MISS_v3$MissSums <- rowSums(MISS_v3 )

#This isolates and sums the Total data count
MISS_v4 <- MISS[,c("Total", "Total.1", "Total.2", "Total.3", "Total.4", "Total.5", "Total.6", "Total.7", "Total.8", "Total.9", "Total.10", "Total.11", "Total.12", "Total.13", "Total.14")]
MISS_v4$TotalSums <- rowSums(MISS_v4)

#This brings together the summed missing and total columns and calculates the proportion of missing data per sample
DF <- as.data.frame(cbind(MISS_v2, MISS_v3$MissSums, MISS_v4$TotalSums))
colnames(DF) <- c("Name", "Missing", "Total")
DF$Missing <- as.numeric(DF$Missing)
DF$Total <- as.numeric(DF$Total)
DF$Perc_Missing <- DF$Missing/DF$Total

#This writes out the final file
write.csv(DF, file="Tassel_TagTaxaDist_AllData_v7_v_K10L2.Sub1Perc.MissingAll.txt", row.names = FALSE, quote = FALSE)

#############Merging all of these separately makes trouble shooting easier.

DF$Simp_Name <- NA
DF$Extra_Name <- NA
DF$OG_Name <- NA

#This section isolates the Swarts Data and merges it with the groups.
GROUPS_SW <- subset(GROUPS, Data_Source == "Swarts_etal_2017")
DF_SW <-merge(DF, GROUPS_SW, by = "Name")

#This section isolates the Romay Data and merges it with the groups.
GROUPS_RY <- subset(GROUPS, Data_Source == "Romay_etal_2013")
DF_RY <- merge(DF, GROUPS_RY, by = "Name")

#This section isolates the Romay Data and merges it with the groups.
GROUPS_DC <- subset(GROUPS, Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2")
DF_DC <- merge(DF, GROUPS_DC, by = "Name")

#This section isolates the Romero Navarro Data and merges it with the groups.
GROUPS_RN <- subset(GROUPS, Data_Source == "Romero-Navarro_etal_2017")
DF_RN_temp <- DF[(grep("RN", DF$Name)),]

#This function splits the Romero Navarro names so they are compatible with the groups file
SPLIT <- function(x){
  DF_RN_temp[x,5] <<- paste(strsplit(DF_RN_temp[x,1], "[.]")[[1]][1], "RN", sep = ".")
  DF_RN_temp[x,6] <<- paste((strsplit(DF_RN_temp[x,1], "[.]")[[1]][2]), (strsplit(DF_RN_temp[x,1], "[.]")[[1]][3]), sep=".")
  DF_RN_temp[x,7] <<- DF_RN_temp[x,1]
}

lapply(c(1:nrow(DF_RN_temp)), SPLIT)

DF_RN <- merge(DF_RN_temp, GROUPS_RN, by.x ="Simp_Name", by.y = "Name")

#This creates a list of all the unique names from the Romero-Navarro dataset
SAMPS <- unique(DF_RN$Simp_Name)

#This pulls just the header as a dataframe
DF_RN_FIX <- DF_RN[1,]
DF_RN_FIX <- DF_RN_FIX[-c(1),]

###################This section sums the missing and total rows for the Romero-Navarro data and recalculates their Proportion Missing
i=1
for(i in 1:length(SAMPS)) {
  x <- SAMPS[i]
  index <- which(DF_RN$Simp_Name %in% x)
  TEMP1 <- as.data.frame(DF_RN[index,3], drop = FALSE)
  SUM_MISS <- colSums(TEMP1)
  TEMP2 <- as.data.frame(DF_RN[index,4], drop = FALSE)
  SUM_TOT <- colSums(TEMP2)
  NEW <- DF_RN[index[1],]
  NEW[1,2] <- x
  NEW[1,3] <- SUM_MISS
  NEW[1,4] <- SUM_TOT
  NEW[1,5] <- SUM_MISS/SUM_TOT
  DF_RN_FIX <<- rbind(DF_RN_FIX, NEW)
}


#This merges all of the data 
DF_GROUPS <- rbind(DF_RN_FIX, DF_SW, DF_RY, DF_DC)

#This plots the distribution of the proportion of missing data by data source
png(filename = "K10L2_MissingDataByDataSource.png", height = 700, width = 500)
ggplot(data=DF_GROUPS) +
  geom_histogram(aes(x=Perc_Missing, fill=Data_Source), bins = 100) +
  facet_wrap(~Data_Source, scales="free_y", ncol = 1) +
  xlab("Proportion of Missing Data By Data Source") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text = element_text(size = 15), legend.position = "none")
dev.off()

DF_GROUPS_v2 <- subset(DF_GROUPS, DF_GROUPS$Data_Source == "Dawe_Lab_1" | DF_GROUPS$Data_Source == "Dawe_Lab_2")

png(filename = "K10L2_MissingDataByAb10Status.png", height = 1100, width = 500)
ggplot(data=DF_GROUPS_v2) +
  geom_histogram(aes(x=Perc_Missing, fill=Ab10_Status), bins = 100) +
  facet_wrap(~Ab10_Status, scales="fixed", ncol = 1) +
  xlab("Proportion of Missing Data By Ab10 Status") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text = element_text(size = 15), legend.position = "none")
dev.off()

#I am determining the lowest threshold I can reasonably keep by taking the lowest proportion of missing data for blank samples and subtracting 0.001. I selected this because it seems reasonable
DF_GROUPS_BLANK <- DF_GROUPS[grep("BLANK", DF_GROUPS$Name),]
CUTOFF <- min(DF_GROUPS_BLANK$Perc_Missing)-0.001

#This applies the cutoff
DF_GROUPS_v3 <- subset(DF_GROUPS, Perc_Missing <= CUTOFF)

#This converts it back to the original separated out Romero-Navarro Data
DF_GROUPS_v3_OTH <- subset(DF_GROUPS_v3, Data_Source != "Romero-Navarro_etal_2017")
DF_GROUPS_v3_RN <- subset(DF_GROUPS_v3, Data_Source == "Romero-Navarro_etal_2017")

UMSUM <- function(x){
  IN <- which(DF_RN$Simp_Name %in% x)
  NEW <- DF_RN[IN,]
  DF_GROUPS_v3_OTH <<- rbind(DF_GROUPS_v3_OTH, NEW)
}

lapply(DF_GROUPS_v3_RN$Simp_Name, UMSUM)


#This writes out the samples that PASSED the filter with their original names
write.csv(DF_GROUPS_v3_OTH, file="K10L2_SamplesPASSEDMissingDataFILTER.csv", row.names = FALSE, quote=FALSE)

#This plots the data after the cutoff
png(filename = "K10L2_FilteredMissingDataByDataSource.png", height = 700, width = 600)
ggplot(data=DF_GROUPS_v3) +
  geom_histogram(aes(x=Perc_Missing, fill=Data_Source), bins = 100) +
  facet_wrap(~Data_Source, scales="free_y", ncol = 1) +
  xlab("Proportion of Missing Data By Data Source") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text = element_text(size = 15), legend.position = "none")
dev.off()

#This isolates and writes out all of the samples that are REMOVED due to missing data
DF_GROUPS_v4 <- subset(DF_GROUPS, Perc_Missing >= CUTOFF)
write.csv(DF_GROUPS_v4, file="K10L2_SamplesRemovedForExcessMissingData.csv", row.names = FALSE, quote = FALSE)

DF_GROUPS_v5 <- subset(DF_GROUPS_v4, Maize_Type != "Blank")

#69 samples that were not blanks were dropped for low coverage
