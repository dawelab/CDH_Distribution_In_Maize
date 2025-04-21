library(glmnet)
library(vcfR)
library(stringr)
library(DHARMa)

setwd("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map")

#Load the data with all of the CDH calls
GROUPS <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumCompleteWholeGenomePCBChr.csv")
colnames(GROUPS)

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

#This drops any variable with a missing BChrom value
DATA <- subset(DATA, is.na(DATA$KMeans_BChrom) == FALSE)

#This drops everything but the BChrom yes no and the variables
SUB <- DATA[,-c(1:14, 16:17, 19:20, 31:115, 135:141)]
#Set the rownames
rownames(SUB) <- DATA$Name
#This drops any row with missing data
SUB_clean <- na.omit(SUB)

#This tests for a relationship to elevation
model1 <- glm(KMeans_BChrom ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + elev, family = "binomial", data = SUB_clean)
summary(model1)
#No assicoation between B chrom and elevation

#I am testing the assumptions
model1_simulationOutput <- simulateResiduals(fittedModel = model1, plot = F)

png("BChr_Elevation_DHARMa_Plots_Start.png")
plot(model1_simulationOutput)
dev.off()
#There are no issues here and the residual deviance is lower than the residual df

#This uses the variables selected via LASSO to run a traditional logisitic regression
model2 <- glm(KMeans_BChrom ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + KMeans_Ab10 + srad + wind + vapr +  sq3_simp + sq4_simp + bio_10 + bio_18 + bio_4 + bio_15 + sq1_simp + sq5_simp + sq6_simp + elev, family = "binomial", data = SUB_clean)

#This calculates the percent of deviance explained
(model2$null.deviance - model2$deviance) / model2$null.deviance * 100

#This summarizes the model
summary(model2)

#I am testing the assumptions
model2_simulationOutput <- simulateResiduals(fittedModel = model2, plot = F)

png("BChr_All_DHARMa_Plots_Start.png")
plot(model2_simulationOutput)
dev.off()
#There are no issues here and the residual deviance is lower than the residual df

#This simplifies the model
step(model2)
#This is robust to rearrangment 

#I began with the model recommended by step and simplified further.
model3 <- glm(KMeans_BChrom ~ srad + sq5_simp + PC1 + PC2 + PC5 + PC8 + PC9 + PC10, family = "binomial", data = SUB_clean)
summary(model3)
#This is robust to rearrangment

#This calculates the percent of deviance explained
(model3$null.deviance - model3$deviance) / model3$null.deviance * 100

#This simplifies the model further by collapsing the soil variables
#SUB_clean$sq5_simp <- gsub("1", "Good", SUB_clean$sq5_simp)
#SUB_clean$sq5_simp <- gsub("2", "Good", SUB_clean$sq5_simp)
#SUB_clean$sq5_simp <- gsub("3", "Bad", SUB_clean$sq5_simp)
#SUB_clean$sq5_simp <- gsub("4", "Very Bad", SUB_clean$sq5_simp)
#SUB_clean$sq5_simp <- factor(SUB_clean$sq5_simp, levels = c("Good", "Bad", "Very Bad"))
#summary(as.factor(SUB_clean$sq5_simp))


#This drops some of the variables that seem non significant 
model4 <- glm(KMeans_BChrom ~ srad + sq5_simp + PC1 + PC2 + PC5 + PC8 + PC9 + PC10, family = "binomial", data = SUB_clean)
#This summarizes the model
summary(model4)
#This calculates the percent of deviance explained
(model4$null.deviance - model4$deviance) / model4$null.deviance * 100

#This tests the assumptions
model4_simulationOutput <- simulateResiduals(fittedModel = model4, plot = F)

png("BChr_All_DHARMa_Plots_End.png")
plot(model4_simulationOutput)
dev.off()


###################################################################################
#This saves the model coefficients for plotting
###################################################################################

#This saves the model coefficients
#This extracts the variable names
Name <- names(coef(model4))[-c(1)]
#This extracts the beta values
Beta <- coef(model4)[-c(1)]
#This extracts the p values
p <- coef(summary(model4))[,4][-c(1)]

SUM <- data.frame("Variable"=Name, "Beta"=Beta, "p"=p)
SUM$CDH <- "BChr"

write.csv(SUM, "BChr_GLM_FinalModel.csv")

###################################################################################
#This attempts to add the associated loci
###################################################################################

#This loads the association data from the GWAS
BChr_assoc <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChrom_results_PCA.assoc.logistic.filt")
BChr_assoc$neglog <- -log10(BChr_assoc$P)

#This filters the variants to only those that are HIGHLY significant
SIG <- subset(BChr_assoc, neglog >= 7.3)

write.table(SIG, "BChr_SigAssocSNP.txt", quote = FALSE, row.names = FALSE)


VCF <- vcfR::read.vcfR("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/GWAS_Mo17_Bchr.chr.landrace.filt.imputed.filt.vcf.gz")

#Exrtact SNP position
FIX <- as.data.frame(VCF@fix, stringsAsFactors = FALSE)

# Extract genotype data
GT <- as.data.frame(VCF@gt, stringsAsFactors = FALSE)

DF <- cbind(FIX, GT)

GE <- DF[which(DF$ID %in% SIG$SNP),]

convert <- function(x) {
  y <- str_split((x), ":")[[1]][1]
  y <- sub("1\\|0", "1", y)
  y <- sub("0\\|1", "1", y)
  y <- sub("1\\|1", "2", y)
  y <- sub("0\\|0", "0", y)
  y <- sub("1\\|2", "3", y)
  y <- sub("2\\|1", "3", y)
  y <- sub("0\\|2", "4", y)
  y <- sub("2\\|0", "4", y)
  y <- sub("2\\|2", "5", y)
}


GE[,-c(1:9)] <- apply(GE[,-c(1:9)], c(1, 2), convert)

# Replace "./." with NA in genotype data
GE[,-c(1:9)] <- lapply(GE[,-c(1:9)], function(x) gsub("\\./\\.", NA, x))


#This reshapes the data to be merged with the rest of the environmental data
rownames(GE) <- GE$ID
GE_SUB <- GE[,-c(1:9)]
GE_t <- as.data.frame(t(GE_SUB))

#This creates columns to merge by
GE_t$Name <- rownames(GE_t)
SUB_clean$Name <- rownames(SUB_clean)
ALL <- merge(SUB_clean, GE_t, by = "Name")

#This drops the Name column used to merge
ALL <- ALL[,-c(1)]

#this sets 0 as a factor level for SNPs
ALL$S3_211286197 <- factor(ALL$S3_211286197, levels=c(0,1,2,3,4,5))
ALL$S4_248470517 <- factor(ALL$S4_248470517, levels=c(0,1,2,3,4,5))


#This drops some of the variables that seem non significant 
model6 <- glm(KMeans_BChrom ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + KMeans_Ab10 + srad + wind + vapr +  sq3_simp + sq4_simp + bio_10 + bio_18 + bio_4 + bio_15 + sq1_simp + sq5_simp + sq6_simp + elev + S3_211286197 + S4_248470517, family = "binomial", data = ALL)
summary(model6)
#This calculates the percent of deviance explained
(model6$null.deviance - model6$deviance) / model6$null.deviance * 100


###################################################################################
#The above demonstrates that none of the 2nd alleles are significantly associated I am removing them
###################################################################################

GE <- DF[which(DF$ID %in% SIG$SNP),]

convert <- function(x) {
  y <- str_split((x), ":")[[1]][1]
  y <- sub("1\\|0", "1", y)
  y <- sub("0\\|1", "1", y)
  y <- sub("1\\|1", "2", y)
  y <- sub("0\\|0", "0", y)
  y <- sub("1\\|2", NA, y)
  y <- sub("2\\|1", NA, y)
  y <- sub("0\\|2", NA, y)
  y <- sub("2\\|0", NA, y)
  y <- sub("2\\|2", NA, y)
}

GE[,-c(1:9)] <- apply(GE[,-c(1:9)], c(1, 2), convert)

# Replace "./." with NA in genotype data
GE[,-c(1:9)] <- lapply(GE[,-c(1:9)], function(x) gsub("\\./\\.", NA, x))


#This reshapes the data to be merged with the rest of the environmental data
rownames(GE) <- GE$ID
GE_SUB <- GE[,-c(1:9)]
GE_t <- as.data.frame(t(GE_SUB))

#This creates columns to merge by
GE_t$Name <- rownames(GE_t)
SUB_clean$Name <- rownames(SUB_clean)
ALL <- merge(SUB_clean, GE_t, by = "Name")

#This writes out the All data frame in case I need it later
write.csv(ALL, "~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChr_All_Variables.csv", row.names = FALSE, quote = FALSE)

ALL <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChr_All_Variables.csv")

#This drops the Name column used to merge
ALL <- ALL[,-c(1)]

#This reclassify the soil variables as factors
ALL$sq1_simp <- as.factor(ALL$sq1_simp)
ALL$sq2_simp <- as.factor(ALL$sq2_simp)
ALL$sq3_simp <- as.factor(ALL$sq3_simp)
ALL$sq4_simp <- as.factor(ALL$sq4_simp)
ALL$sq5_simp <- as.factor(ALL$sq5_simp)
ALL$sq6_simp <- as.factor(ALL$sq6_simp)


#this sets 0 as a factor level for SNPs
ALL$S3_211286197 <- as.numeric(ALL$S3_211286197)
ALL$S4_248470517 <- as.numeric(ALL$S4_248470517)

#This selects only complete cases
ALL <- ALL[complete.cases(ALL), ]

#This drops some of the variables that seem non significant 
model7 <- glm(KMeans_BChrom ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + KMeans_Ab10 + srad + wind + vapr +  sq3_simp + sq4_simp + bio_10 + bio_18 + bio_4 + bio_15 + sq1_simp + sq5_simp + sq6_simp + S3_211286197 + S4_248470517 + elev, family = "binomial", data = ALL)
summary(model7)
#This calculates the percent of deviance explained
(model7$null.deviance - model7$deviance) / model7$null.deviance * 100

model7_simulationOutput <- simulateResiduals(fittedModel = model7, plot = F)

png("BChr_AllSNPs_DHARMa_Plots_Start.png")
plot(model7_simulationOutput)
dev.off()

step(model7)
#This is robust to rearrangement 


#I began with the model recommended by step and simplified further
model8 <- glm(KMeans_BChrom ~ PC1 + PC2 + PC5 + PC8 + PC9 + PC10 + 
                srad + sq5_simp + S3_211286197 + 
                S4_248470517, family = "binomial", data = ALL)
summary(model8)
#This is robust to rearrangement 

#This simplifies the model further by collapsing the soil variables
ALL$sq5_simp <- gsub("1", "Good", ALL$sq5_simp)
ALL$sq5_simp <- gsub("2", "Good", ALL$sq5_simp)
ALL$sq5_simp <- gsub("3", "Bad", ALL$sq5_simp)
ALL$sq5_simp <- gsub("4", "Bad", ALL$sq5_simp)
ALL$sq5_simp <- factor(ALL$sq5_simp, levels = c("Good", "Bad"))
summary(as.factor(ALL$sq5_simp))

#I began with the model recommended by step and simplified further
model8 <- glm(KMeans_BChrom ~ S3_211286197 + S4_248470517 + PC1 + PC2 + PC5 + PC8 + PC9 + PC10 + srad, family = "binomial", data = ALL)
summary(model8)
#This is robust to rearrangement 


#This calculates the percent of deviance explained
(model8$null.deviance - model8$deviance) / model8$null.deviance * 100

#This tests assumptions
model8_simulationOutput <- simulateResiduals(fittedModel = model8, plot = F)

png("BChr_AllSNPs_DHARMa_Plots_End.png")
plot(model8_simulationOutput)
dev.off()
#There are no issues here and the residual deviance is lower than the residual df

###################################################################################
#This saves the model coefficients for plotting
###################################################################################

#This extracts the variable names
Name <- names(coef(model8))[-c(1)]
#This extracts the beta values
Beta <- coef(model8)[-c(1)]
#This extracts the p values
p <- coef(summary(model8))[,4][-c(1)]

SUM <- data.frame("Variable"=Name, "Beta"=Beta, "p"=p)
SUM$CDH <- "BChr"

write.csv(SUM, "BChr_GLM_FinalModelSNP.csv")


###################################################################################
#This does variance partitioning
###################################################################################

#In order for this to work I need to drop anything with missing data. I loose over half the data this way
ALL_Comp <- na.omit(ALL) 

#Several Variables drop out I removed them stepwise based on significance 
modelFULL <- glm(KMeans_BChrom ~ PC1 + PC2 + PC5 + PC8 + PC9 + PC10 + srad + S3_211286197 + S4_248470517, family = "binomial", data = ALL_Comp)
summary(modelFULL)
#This calculates the percent of deviance explained
(modelFULL$null.deviance - modelFULL$deviance) / modelFULL$null.deviance * 100


modelPC <- glm(KMeans_BChrom ~ PC1 + PC2 + PC5 + PC8 + PC9 + PC10, family = "binomial", data = ALL_Comp)
summary(modelPC)
#This calculates the percent of deviance explained
(modelPC$null.deviance - modelPC$deviance) / modelPC$null.deviance * 100


modelEV <- glm(KMeans_BChrom ~ srad, family = "binomial", data = ALL_Comp)
summary(modelEV)
#This calculates the percent of deviance explained
(modelEV$null.deviance - modelEV$deviance) / modelEV$null.deviance * 100


modelSNP <- glm(KMeans_BChrom ~ S3_211286197 + S4_248470517, family = "binomial", data = ALL_Comp)
summary(modelSNP)
#This calculates the percent of deviance explained
(modelSNP$null.deviance - modelSNP$deviance) / modelSNP$null.deviance * 100


#This creates a report on the variance partitioning
Report <- data.frame("CDH" = "B Chr","Full_Model"= (modelFULL$null.deviance - modelFULL$deviance) / modelFULL$null.deviance * 100, "Pop_Structure"= (modelPC$null.deviance - modelPC$deviance) / modelPC$null.deviance * 100, "Env"= (modelEV$null.deviance - modelEV$deviance) / modelEV$null.deviance * 100, "SNPs"=(modelSNP$null.deviance - modelSNP$deviance) / modelSNP$null.deviance * 100)

#This reads in the report
Report_Full <- read.csv("CDH_DeviancePartitioning.csv")

#This binds the Ab10 sample with rest of the report
Report_Full <- rbind(Report_Full, Report)

#This writes the report
write.csv(Report_Full, "CDH_DeviancePartitioning.csv", row.names = FALSE, quote = FALSE)
