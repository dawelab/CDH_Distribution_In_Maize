
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
DATA$KMeans_Ab10 <- gsub("Ambiguous", NA, DATA$KMeans_Ab10)
DATA$KMeans_Ab10 <- gsub("Ab10", 1, DATA$KMeans_Ab10)
DATA$KMeans_Ab10 <- gsub("N10", 0, DATA$KMeans_Ab10)
DATA$KMeans_Ab10 <- as.factor(DATA$KMeans_Ab10)
DATA$sq1 <- factor(DATA$sq1, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq2 <- factor(DATA$sq1, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq3 <- factor(DATA$sq3, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq4 <- factor(DATA$sq4, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq5 <- factor(DATA$sq5, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq6 <- factor(DATA$sq6, levels = c("1", "2", "3", "4", "5", "6", "7"))

DATA$sq1_simp <- gsub("1", "Good", DATA$sq1)
DATA$sq1_simp <- gsub("2", "Good", DATA$sq1)
DATA$sq1_simp <- gsub("3", "Bad", DATA$sq1)
DATA$sq1_simp <- gsub("4", "Bad", DATA$sq1)
DATA$sq1_simp <- gsub("5", NA, DATA$sq1)
DATA$sq1_simp <- gsub("6", NA, DATA$sq1)
DATA$sq1_simp <- gsub("7", NA, DATA$sq1)

DATA$sq2_simp <- gsub("1", "Good", DATA$sq2)
DATA$sq2_simp <- gsub("2", "Good", DATA$sq2)
DATA$sq2_simp <- gsub("3", "Bad", DATA$sq2)
DATA$sq2_simp <- gsub("4", "Bad", DATA$sq2)
DATA$sq2_simp <- gsub("5", NA, DATA$sq2)
DATA$sq2_simp <- gsub("6", NA, DATA$sq2)
DATA$sq2_simp <- gsub("7", NA, DATA$sq2)

DATA$sq3_simp <- gsub("1", "Good", DATA$sq3)
DATA$sq3_simp <- gsub("2", "Good", DATA$sq3)
DATA$sq3_simp <- gsub("3", "Bad", DATA$sq3)
DATA$sq3_simp <- gsub("4", "Bad", DATA$sq3)
DATA$sq3_simp <- gsub("5", NA, DATA$sq3)
DATA$sq3_simp <- gsub("6", NA, DATA$sq3)
DATA$sq3_simp <- gsub("7", NA, DATA$sq3)

DATA$sq4_simp <- gsub("1", "Good", DATA$sq4)
DATA$sq4_simp <- gsub("2", "Good", DATA$sq4)
DATA$sq4_simp <- gsub("3", "Bad", DATA$sq4)
DATA$sq4_simp <- gsub("4", "Bad", DATA$sq4)
DATA$sq4_simp <- gsub("5", NA, DATA$sq4)
DATA$sq4_simp <- gsub("6", NA, DATA$sq4)
DATA$sq4_simp <- gsub("7", NA, DATA$sq4)

DATA$sq5_simp <- gsub("1", "Good", DATA$sq5)
DATA$sq5_simp <- gsub("2", "Good", DATA$sq5)
DATA$sq5_simp <- gsub("3", "Bad", DATA$sq5)
DATA$sq5_simp <- gsub("4", "Bad", DATA$sq5)
DATA$sq5_simp <- gsub("5", NA, DATA$sq5)
DATA$sq5_simp <- gsub("6", NA, DATA$sq5)
DATA$sq5_simp <- gsub("7", NA, DATA$sq5)

DATA$sq6_simp <- as.factor(gsub("1", "Good", DATA$sq6))
DATA$sq6_simp <- gsub("2", "Good", DATA$sq6)
DATA$sq6_simp <- gsub("3", "Bad", DATA$sq6)
DATA$sq6_simp <- gsub("4", "Bad", DATA$sq6)
DATA$sq6_simp <- gsub("5", NA, DATA$sq6)
DATA$sq6_simp <- gsub("6", NA, DATA$sq6)
DATA$sq6_simp <- gsub("7", NA, DATA$sq6)

#This drops any variable with a missing Ab10 value
DATA <- subset(DATA, is.na(DATA$KMeans_Ab10) == FALSE)

#This tests for a relationship between Ab10 and elevation
model1 <- glm(KMeans_Ab10 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + elev, family = "binomial", data=DATA)
summary(model1)

#This tests for a relationship between Ab10 and and any of the non co linear environmental variables
model2 <- glm(KMeans_Ab10 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + srad_03 + srad_07 + srad_10 + wind_06 + bio_10 + bio_12 + bio_15 + bio_18 + bio_19 + bio_2 + bio_4 + sq1_simp + sq3_simp + sq4_simp + sq5_simp + sq6_simp, family = "binomial", data=DATA)
summary(model2)

model2 <- glm(KMeans_Ab10 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + + srad_03 + srad_07 + srad_10 + wind_06 + bio_1 + bio_3 + bio_1 + bio_2 + + bio_4 + + bio_4 + + bio_5 + + bio_6 + + bio_7 + bio_8 + bio_9 + bio_10 + bio_11 + bio_12 + bio_13 + bio_14 + bio_15 + bio_16 + bio_17 + bio_18 + bio_19 + sq1_simp + sq2_simp +sq3_simp+ sq4_simp +sq5_simp+ sq6_simp, family = "binomial", data=DATA)
summary(model2)
#This calculates the percent of deiance explained
(model2$null.deviance - model2$deviance) / model2$null.deviance * 100

#This automaticallt attempts to simplify the model 
step(model2, test="Chi")

model_test <- glm(formula = KMeans_Ab10 ~ PC1 + PC2 + PC10 + 
                    srad_03 + srad_07 + bio_1 + bio_4 + bio_5 + 
                    bio_6 + bio_7 + bio_9 + bio_11 +   sq1_simp + 
                    sq5_simp + sq6_simp, family = "binomial", data = DATA)

#This calculates the percent of deiance explained
(model_test$null.deviance - model_test$deviance) / model_test$null.deviance * 100
summary(model_test)

model4 <- glm(KMeans_Ab10 ~ PC1 + PC2 + PC3 + PC5 + PC6 + srad_07 + wind_06 + bio_15 + bio_4 + sq1 + sq5 + sq6, family = "binomial", data=DATA)
model4

drop1(model4, test="Chi")

model5 <- glm(KMeans_Ab10 ~ PC1 + PC2 + PC3 + PC5 + PC6 + srad_07 + wind_06 + bio_4 + sq1 + sq5 + sq6, family = "binomial", data=DATA)
model5
drop1(model5, test="Chi")

model6 <- glm(KMeans_Ab10 ~ PC1 + PC2 + PC3 + PC5 + PC6 + wind_06 + bio_4 + sq1 + sq5 + sq6, family = "binomial", data=DATA)
model6
drop1(model6, test="Chi")
summary(model6)
