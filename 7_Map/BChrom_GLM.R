
#Load the data with all of the CDH calls
GROUPS <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumCompleteWholeGenomePC.csv")
colnames(GROUPS)

#Load the environmental data for all the GPS points from WorldClim2 and FAO
ENV <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")
colnames(ENV)

#This brings together the calls and the environmental data
DATA <- merge(GROUPS, ENV, by=c("Name", "Data_Source", "Accession", "Group", "Maize_Type", "Ab10_Status", "B_Chrom_Status", "Latitude", "Longitude", "Altitude", "DTA_BLUP", "DTS_BLUP", "PH_BLUP", "PresTiller_BLUP"))

colnames(DATA) <- gsub("wc2.1_30s_", "", colnames(DATA))

#This filters to only lines with environmental data
DATA <- subset(DATA, is.na(DATA$Latitude) == FALSE )

#This reformats the data
DATA$KMeans_BChrom <- gsub("Ambiguous", NA, DATA$KMeans_BChrom )
DATA$KMeans_BChrom <- gsub("Yes", 1, DATA$KMeans_BChrom )
DATA$KMeans_BChrom <- gsub("No", 0, DATA$KMeans_BChrom )
DATA$KMeans_BChrom <- as.factor(DATA$KMeans_BChrom )
DATA$sq1 <- factor(DATA$sq1, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq3 <- factor(DATA$sq3, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq4 <- factor(DATA$sq4, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq5 <- factor(DATA$sq5, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq6 <- factor(DATA$sq6, levels = c("1", "2", "3", "4", "5", "6", "7"))

#This drops any variable with a missing BChrom value
DATA <- subset(DATA, is.na(DATA$KMeans_BChrom) == FALSE)

#This tests for a relationship between BChrom and elevation
model1 <- glm(KMeans_BChrom  ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +  elev, family = "binomial", data=DATA)

summary(model1)

#This tests for a relationship between BChrom and and any of the non co linear environmental variables
model2 <- glm(KMeans_BChrom ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + srad_03 + srad_07 + srad_10 + wind_06 + bio_10 + bio_12 + bio_15 + bio_18 + bio_19 + bio_2 + bio_4 + sq1 + sq3 + sq4 + sq5 + sq6, family = "binomial", data=DATA)

summary(model2)

#This looks for a relationship between B chromosome pseudo copy number and elevation

#This subsets to only B chromosome containing lines
DATA <- subset(DATA, KMeans_BChrom == 1)

model3 <- glm(Mean  ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + elev, family = "gaussian", data=DATA)

summary(model3)

#This does the same for all of the non co linear environmental variables
model4 <- glm(Mean ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + srad_03 + srad_07 + srad_10 + wind_06 + bio_10 + bio_12 + bio_15 + bio_18 + bio_19 + bio_2 + bio_4 + sq1 + sq3 + sq4 + sq5 + sq6, family = "gaussian", data=DATA)

summary(model4)
