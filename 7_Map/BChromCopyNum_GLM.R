library(glmnet)
library(vcfR)
library(stringr)
library(DHARMa)

#Load the data with all of the CDH calls
GROUPS <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumCompleteWholeGenomePCBChr.csv")
colnames(GROUPS)

#This loads in B chromosome pseudo copy number
Copy <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/6_ClassifyAllData/BChr_ContExperimental_CopyNumber.csv")
Copy <- Copy[,c("Name", "Mean")]
colnames(Copy) <- c("Name", "BChrom_Pseudo_Copy_Number")

#Load the environmental data for all the GPS points from WorldClim2 and FAO
ENV <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")
colnames(ENV)

#This brings together the calls and the environmental data
DATA <- merge(GROUPS, ENV, by=c("Name", "Data_Source", "Accession", "Group", "Maize_Type", "Ab10_Status", "B_Chrom_Status", "Latitude", "Longitude", "Altitude", "DTA_BLUP", "DTS_BLUP", "PH_BLUP", "PresTiller_BLUP"))

colnames(DATA) <- gsub("wc2.1_30s_", "", colnames(DATA))

#This filters to only lines with environmental data
DATA <- subset(DATA, is.na(DATA$Latitude) == FALSE )

###################################################################################
#This reshapes the data
###################################################################################

#This generates PCs from solar radiation, wind, and vapor pressure
srad <- DATA[,c(79:90)]
srad$sradAvg <- rowMeans(srad)
srad <- srad[,13]

wind <- DATA[,c(91:102)]
wind$windAvg <- rowMeans(wind)
wind <- wind[,13]

vapr <- DATA[,c(103:114)]
vapr$vaprAvg <- rowMeans(vapr)
vapr <- vapr[,13]

#This adds the PC columns to the DATA
DATA <- cbind(DATA, srad, wind, vapr)

#This reformats the data
DATA$KMeans_BChrom <- gsub("Ambiguous", NA, DATA$KMeans_BChrom)
DATA$KMeans_BChrom <- gsub("Yes", 1, DATA$KMeans_BChrom)
DATA$KMeans_BChrom <- gsub("No", 0, DATA$KMeans_BChrom)
DATA$KMeans_BChrom <- as.factor(DATA$KMeans_BChrom)

DATA$KMeans_Ab10 <- gsub("Ambiguous", NA, DATA$KMeans_Ab10)
DATA$KMeans_Ab10 <- gsub("Ab10", 1, DATA$KMeans_Ab10)
DATA$KMeans_Ab10 <- gsub("N10", 0, DATA$KMeans_Ab10)
DATA$KMeans_Ab10 <- as.factor(DATA$KMeans_Ab10)

DATA$KMeans_K10L2 <- gsub("Ambiguous", NA, DATA$KMeans_K10L2)
DATA$KMeans_K10L2 <- gsub("K10L2", 1, DATA$KMeans_K10L2)
DATA$KMeans_K10L2 <- gsub("N10", 0, DATA$KMeans_K10L2)
DATA$KMeans_K10L2 <- as.factor(DATA$KMeans_K10L2)


DATA$sq1 <- factor(DATA$sq1, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq2 <- factor(DATA$sq2, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq3 <- factor(DATA$sq3, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq4 <- factor(DATA$sq4, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq5 <- factor(DATA$sq5, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq6 <- factor(DATA$sq6, levels = c("1", "2", "3", "4", "5", "6", "7"))

#DATA$sq1_simp <- gsub("1", "Good", DATA$sq1)
#DATA$sq1_simp <- gsub("2", "Good", DATA$sq1_simp)
#DATA$sq1_simp <- gsub("3", "Bad", DATA$sq1_simp)
#DATA$sq1_simp <- gsub("4", "Bad", DATA$sq1_simp)
DATA$sq1_simp <- gsub("5", NA, DATA$sq1)
DATA$sq1_simp <- gsub("6", NA, DATA$sq1_simp)
DATA$sq1_simp <- gsub("7", NA, DATA$sq1_simp)
DATA$sq1_simp <- factor(DATA$sq1_simp, levels = c("1", "2", "3", "4"))

#DATA$sq2_simp <- gsub("1", "Good", DATA$sq2)
#DATA$sq2_simp <- gsub("2", "Good", DATA$sq2_simp)
#DATA$sq2_simp <- gsub("3", "Bad", DATA$sq2_simp)
#DATA$sq2_simp <- gsub("4", "Bad", DATA$sq2_simp)
DATA$sq2_simp <- gsub("5", NA, DATA$sq2)
DATA$sq2_simp <- gsub("6", NA, DATA$sq2_simp)
DATA$sq2_simp <- gsub("7", NA, DATA$sq2_simp)
DATA$sq2_simp <- factor(DATA$sq2_simp, levels = c("1", "2", "3", "4"))

#DATA$sq3_simp <- gsub("1", "Good", DATA$sq3)
#DATA$sq3_simp <- gsub("2", "Good", DATA$sq3_simp)
#DATA$sq3_simp <- gsub("3", "Bad", DATA$sq3_simp)
#DATA$sq3_simp <- gsub("4", "Bad", DATA$sq3_simp)
DATA$sq3_simp <- gsub("5", NA, DATA$sq3)
DATA$sq3_simp <- gsub("6", NA, DATA$sq3_simp)
DATA$sq3_simp <- gsub("7", NA, DATA$sq3_simp)
DATA$sq3_simp <- factor(DATA$sq3_simp, levels = c("1", "2", "3", "4"))

#DATA$sq4_simp <- gsub("1", "Good", DATA$sq4)
#DATA$sq4_simp <- gsub("2", "Good", DATA$sq4_simp)
#DATA$sq4_simp <- gsub("3", "Bad", DATA$sq4_simp)
#DATA$sq4_simp <- gsub("4", "Bad", DATA$sq4_simp)
DATA$sq4_simp <- gsub("5", NA, DATA$sq4)
DATA$sq4_simp <- gsub("6", NA, DATA$sq4_simp)
DATA$sq4_simp <- gsub("7", NA, DATA$sq4_simp)
DATA$sq4_simp <- factor(DATA$sq4_simp, levels = c("1", "2", "3", "4"))

#DATA$sq5_simp <- gsub("1", "Good", DATA$sq5)
#DATA$sq5_simp <- gsub("2", "Good", DATA$sq5_simp)
#DATA$sq5_simp <- gsub("3", "Bad", DATA$sq5_simp)
#DATA$sq5_simp <- gsub("4", "Bad", DATA$sq5_simp)
DATA$sq5_simp <- gsub("5", NA, DATA$sq5)
DATA$sq5_simp <- gsub("6", NA, DATA$sq5_simp)
DATA$sq5_simp <- gsub("7", NA, DATA$sq5_simp)
DATA$sq5_simp <- factor(DATA$sq5_simp, levels = c("1", "2", "3", "4"))

#DATA$sq6_simp <- gsub("1", "Good", DATA$sq6)
#DATA$sq6_simp <- gsub("2", "Good", DATA$sq6_simp)
#DATA$sq6_simp <- gsub("3", "Bad", DATA$sq6_simp)
#DATA$sq6_simp <- gsub("4", "Bad", DATA$sq6_simp)
DATA$sq6_simp <- gsub("5", NA, DATA$sq6)
DATA$sq6_simp <- gsub("6", NA, DATA$sq6_simp)
DATA$sq6_simp <- gsub("7", NA, DATA$sq6_simp)
DATA$sq6_simp <- factor(DATA$sq6_simp, levels = c("1", "2", "3", "4"))

#This selects only points that have B chromosomes 
DATA <- subset(DATA, DATA$KMeans_BChrom == 1)

#This drops the mean column. It is old and I think wrong
DATA <- DATA[,-c(20)]

#This brings in the updated B chrom copy number
DATA <- merge(DATA, Copy, by = "Name")

#This drops everything but the Ab10 yes no and the variables
SUB <- DATA[,-c(1:14, 16, 19, 31:113, 134:140)]

#Set the rownames
rownames(SUB) <- DATA$Name

#I am not dropping NA at this time becase I want to maintain both Ab10 and K10L2 lines
SUB_clean <- SUB

#I am transforming my data to a normal distribution
SUB_clean$BChrom_Pseudo_Copy_Number <- log10(SUB_clean$BChrom_Pseudo_Copy_Number)

#This tests for a relationship to elevation
model1 <- glm(BChrom_Pseudo_Copy_Number ~ elev + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "gaussian", data = SUB_clean)
summary(model1)
#No relationship between elevation and B chromosome copy number

#I am testing the assumptions
model1_simulationOutput <- simulateResiduals(fittedModel = model1, plot = F)

png("BopyBChr_Elevation_DHARMa_Plots.png")
plot(model1_simulationOutput)
dev.off()
#There are only mild issues here and the residual deviance is less than the degrees of freedom


#This is the largest initial model 
model2 <- glm(BChrom_Pseudo_Copy_Number ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + srad + KMeans_Ab10 + wind + vapr +  sq3_simp + sq4_simp + elev + bio_10 + bio_18 + bio_4 + bio_15 + sq1_simp + sq5_simp + sq6_simp, family = "gaussian", data = SUB_clean)

#This summarizes the model
summary(model2)
#This calculates the percent of deviance explained
(model2$null.deviance - model2$deviance) / model2$null.deviance * 100

#I am testing the assumptions
model2_simulationOutput <- simulateResiduals(fittedModel = model2, plot = F)

png("BopyBChr_All_DHARMa_Plots_Start.png")
plot(model1_simulationOutput)
dev.off()

#There is no relationship between Ab10 and B chromosome copy number. This removes any Ab10 lines so that I can test for K10L2
SUB_clean_K10L2 <- subset(SUB_clean, KMeans_Ab10 == 0)

model2 <- glm(BChrom_Pseudo_Copy_Number ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + srad + KMeans_K10L2 + wind + vapr +  sq3_simp + sq4_simp + elev + bio_10 + bio_18 + bio_4 + bio_15 + sq1_simp + sq5_simp + sq6_simp, family = "gaussian", data = SUB_clean_K10L2)

#This summarizes the model
summary(model2)
#This calculates the percent of deviance explained
(model2$null.deviance - model2$deviance) / model2$null.deviance * 100

#There is no relationship between K10L2 and B chromosome copy number

#I am dropping K10L2 and removing missing data
SUB_clean <- na.omit(SUB_clean)

#This excludes Ab10 and K10L2 for simplicity
model2 <- glm(BChrom_Pseudo_Copy_Number ~ srad + wind + vapr +  sq3_simp + sq4_simp + elev + bio_10 + bio_18 + bio_4 + bio_15 + sq1_simp + sq5_simp + sq6_simp + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "gaussian", data = SUB_clean)


#This simplifies the model
step(model2)
#This is robust to rearangment

#I began with the model recommended by step and simplified from there
model3 <- glm(BChrom_Pseudo_Copy_Number ~ PC4 + vapr + bio_10 + bio_4, family = "gaussian", data = SUB_clean)
summary(model3)
#This is robust to rearrangment
#This calculates the percent of deviance explained
(model3$null.deviance - model3$deviance) / model3$null.deviance * 100

#I am testing the assumptions
model3_simulationOutput <- simulateResiduals(fittedModel = model3, plot = F)

png("BopyBChr_All_DHARMa_Plots_End.png")
plot(model3_simulationOutput)
dev.off()
#There are only slight issues here and the residual deviance is lower than the df

###################################################################################
#This saves the model coefficients for plotting
###################################################################################

#This saves the model coefficients
#This extracts the variable names
Name <- names(coef(model3))[-c(1)]
#This extracts the beta values
Beta <- coef(model3)[-c(1)]
#This extracts the p values
p <- coef(summary(model3))[,4][-c(1)]

SUM <- data.frame("Variable"=Name, "Beta"=Beta, "p"=p)
SUM$CDH <- "BChr Copy Numer"

write.csv(SUM, "BChrCopyNum_GLM_FinalModel.csv")


###################################################################################
#This does variance partitioning
###################################################################################

#In order for this to work I need to drop anything with missing data. I loose over half the data this way
ALL_Comp <- na.omit(SUB_clean) 

#Several Variables drop out I removed them stepwise based on significance 
modelFULL <- glm(BChrom_Pseudo_Copy_Number ~ PC4 + vapr + bio_10 + bio_4, family = "gaussian", data = ALL_Comp)
summary(modelFULL)
#This calculates the percent of deviance explained
(modelFULL$null.deviance - modelFULL$deviance) / modelFULL$null.deviance * 100

modelPC <- glm(BChrom_Pseudo_Copy_Number ~ PC4, family = "gaussian", data = ALL_Comp)
summary(modelPC)
#This calculates the percent of deviance explained
(modelPC$null.deviance - modelPC$deviance) / modelPC$null.deviance * 100

modelEV <- glm(BChrom_Pseudo_Copy_Number ~ bio_10 + bio_4, family = "gaussian", data = ALL_Comp)
summary(modelEV)
#This calculates the percent of deviance explained
(modelEV$null.deviance - modelEV$deviance) / modelEV$null.deviance * 100

#This creates a report on the variance partitioning
Report <- data.frame("CDH" = "B Chr Copy Number","Full_Model"= (modelFULL$null.deviance - modelFULL$deviance) / modelFULL$null.deviance * 100, "Pop_Structure"= (modelPC$null.deviance - modelPC$deviance) / modelPC$null.deviance * 100, "Env"= (modelEV$null.deviance - modelEV$deviance) / modelEV$null.deviance * 100, "SNPs"=NA)

#This reads in the report
Report_Full <- read.csv("CDH_DeviancePartitioning.csv")

#This binds the Ab10 sample with rest of the report
Report_Full <- rbind(Report_Full, Report)

#This writes out the full report
write.csv(Report_Full, "CDH_DeviancePartitioning.csv", row.names = FALSE, quote = FALSE)
