library(glmnet)
library(vcfR)
library(stringr)
library(DHARMa)

setwd("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map")

#Load the data with all of the CDH calls
GROUPS <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumCompleteWholeGenomePCK10L2.csv")
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
DATA$KMeans_K10L2 <- gsub("Ambiguous", NA, DATA$KMeans_K10L2)
DATA$KMeans_K10L2 <- gsub("K10L2", 1, DATA$KMeans_K10L2)
DATA$KMeans_K10L2 <- gsub("N10", 0, DATA$KMeans_K10L2)
DATA$KMeans_K10L2 <- as.factor(DATA$KMeans_K10L2)

DATA$KMeans_BChrom <- gsub("Ambiguous", NA, DATA$KMeans_BChrom)
DATA$KMeans_BChrom <- gsub("Yes", 1, DATA$KMeans_BChrom)
DATA$KMeans_BChrom <- gsub("No", 0, DATA$KMeans_BChrom)
DATA$KMeans_BChrom <- as.factor(DATA$KMeans_BChrom)


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

###################################################################################
#This models the relationship between the CDH and the environment
###################################################################################

#This drops any variable with a missing K10L2 value
DATA <- subset(DATA, is.na(DATA$KMeans_K10L2) == FALSE)

#This drops everything but the K10L2 yes no and the variables
SUB <- DATA[,-c(1:16, 19:20, 31:115, 135:141)]
#Set the rownames
rownames(SUB) <- DATA$Name
#This drops any row with missing data
SUB_clean <- na.omit(SUB)

#This tests for a relationship between K10L2 and elevation
model1 <- glm(KMeans_K10L2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +  elev, family = "binomial", data=DATA)
summary(model1)
#This calculates the percent of deiance explained
(model1$null.deviance - model1$deviance) / model1$null.deviance * 100

#I am testing the assumptions
model1_simulationOutput <- simulateResiduals(fittedModel = model1, plot = F)

png("K10L2_Elevation_DHARMa_Plots.png")
plot(model1_simulationOutput)
dev.off()
#There are no issues here and the null deviance is lower than the degrees of freedom

#This uses the to run a traditional logistic regression
model2 <- glm(KMeans_K10L2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + wind + vapr +  sq3_simp + sq4_simp + elev + bio_10 + bio_18 + bio_4 + bio_15 + sq1_simp + sq5_simp + sq6_simp + KMeans_BChrom + srad, family = "binomial", data = SUB_clean)

#This summarizes the model
summary(model2)

#I am testing the assumptions
model2_simulationOutput <- simulateResiduals(fittedModel = model2, plot = F)

png("K10L2_All_DHARMa_Plots_Start.png")
plot(model2_simulationOutput)
dev.off()
#There are no issues here residual deviance is lower than the residual df


#This calculates the percent of deviance explained
(model2$null.deviance - model2$deviance) / model2$null.deviance * 100

#This simplifies the model based on model AIC
step(model2)
#This model is robust to order

#This drops some of the variables that seem non significant 
model3 <- glm(KMeans_K10L2 ~ vapr + sq3_simp + sq4_simp + bio_10 + PC2 + PC3 + PC10, family = "binomial", data = SUB_clean)
summary(model3)
#This is robust to rearrangement 

#SUB_clean$sq3_simp <- gsub("1", "Good", SUB_clean$sq3_simp)
#SUB_clean$sq3_simp <- gsub("2", "OK", SUB_clean$sq3_simp)
#SUB_clean$sq3_simp <- gsub("3", "Bad", SUB_clean$sq3_simp)
#SUB_clean$sq3_simp <- gsub("4", "Bad", SUB_clean$sq3_simp)
#SUB_clean$sq3_simp <- factor(SUB_clean$sq3_simp, levels = c("Good", "OK", "Bad"))
#summary(as.factor(SUB_clean$sq3_simp))

#SUB_clean$sq4_simp <- gsub("1", "Good", SUB_clean$sq4_simp)
#SUB_clean$sq4_simp <- gsub("2", "Good", SUB_clean$sq4_simp)
#SUB_clean$sq4_simp <- gsub("3", "Bad", SUB_clean$sq4_simp)
#SUB_clean$sq4_simp <- gsub("4", "Bad", SUB_clean$sq4_simp)
#SUB_clean$sq4_simp <- factor(SUB_clean$sq4_simp, levels = c("Good", "Bad"))
#summary(as.factor(SUB_clean$sq4_simp))

#This drops some of the variables that seem non significant 
model4 <- glm(KMeans_K10L2 ~ vapr + sq3_simp + sq4_simp + bio_10 + PC2 + PC3 + PC10, family = "binomial", data = SUB_clean)
summary(model4)
#This calculates the percent of deiance explained
(model4$null.deviance - model4$deviance) / model4$null.deviance * 100

#I am testing the assumptions
model4_simulationOutput <- simulateResiduals(fittedModel = model4, plot = F)

png("K10L2_All_DHARMa_Plots_End.png")
plot(model4_simulationOutput)
dev.off()
#There is some minor deviance here, I'm not worried about it. The residual deviance is lower than the residual degrees of freedom


###################################################################################
#This saves the model coefficients for plotting
###################################################################################

#This extracts the variable names
Name <- names(coef(model4))[-c(1)]
#This extracts the beta values
Beta <- coef(model4)[-c(1)]
#This extracts the p values
p <- coef(summary(model4))[,4][-c(1)]

SUM <- data.frame("Variable"=Name, "Beta"=Beta, "p"=p)
SUM$CDH <- "K10L2"

write.csv(SUM, "K10L2_GLM_FinalModel.csv")

###################################################################################
#This attempts to add the associated loci
###################################################################################
#This loads the association data from the GWAS
K10L2_assoc <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/K10L2_results_PCA.assoc.logistic.filt")
K10L2_assoc$neglog <- -log10(K10L2_assoc$P)

#This filters the variants to only those that are HIGHLY significant
SIG <- subset(K10L2_assoc, neglog >= 7.3)

write.table(SIG, "K10L2_SigAssocSNP.txt", quote = FALSE, row.names = FALSE)


VCF <- vcfR::read.vcfR("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/GWAS_Mo17_K10L2.chr.landrace.filt.imputed.vcf.gz")

#Exrtact SNP position
FIX <- as.data.frame(VCF@fix, stringsAsFactors = FALSE)

# Extract genotype data
GT <- as.data.frame(VCF@gt, stringsAsFactors = FALSE)

DF <- cbind(FIX, GT)

GE <- DF[which(DF$ID %in% SIG$SNP),]

#I found that K10L2 is associated only with 1 and 2 when run as a factor so I'm dropping the second allele and converting them to a numerical

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

rownames(GE) <- GE$ID
GE_SUB <- GE[,-c(1:9)]
GE_t <- as.data.frame(t(GE_SUB))

#This creates columns to merge by
GE_t$Name <- rownames(GE_t)
SUB_clean$Name <- rownames(SUB_clean)
ALL <- merge(SUB_clean, GE_t, by = "Name")

#This converts the one SNP to a numeric
ALL$S1_277994336 <- factor(ALL$S1_277994336, levels=c(0,1,2,3,4,5))
ALL$S4_2010582 <- factor(ALL$S4_2010582, levels=c(0,1,2,3,4,5))
ALL$S4_5674254 <- factor(ALL$S4_5674254, levels=c(0,1,2,3,4,5))
ALL$S5_190520703 <- factor(ALL$S5_190520703 , levels=c(0,1,2,3,4,5))
ALL$S6_175749381 <- factor(ALL$S6_175749381, levels=c(0,1,2,3,4,5))
ALL$S8_163979957 <- factor(ALL$S8_163979957, levels=c(0,1,2,3,4,5))
ALL$S8_164038154 <- factor(ALL$S8_164038154, levels=c(0,1,2,3,4,5))

#This drops the Name column used to merge
ALL <- ALL[,-c(1)]

#This drops some of the variables that seem non significant 
model5 <- glm(KMeans_K10L2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + wind + vapr +  sq3_simp + sq4_simp + elev + bio_10 + bio_18 + bio_4 + bio_15 + sq1_simp + sq5_simp + sq6_simp + KMeans_BChrom + srad + S1_277994336 + S4_2010582 + S4_5674254 + S5_190520703 + S6_175749381 + S8_163979957 + S8_164038154, family = "binomial", data = ALL)
summary(model5)

#This calculates the percent of deviance explained
(model5$null.deviance - model5$deviance) / model5$null.deviance * 100
#This accounts for 14.93% of deviance


##################################################################################################################
#This revealed that the second alleles present were not significantly associated with Ab10 I am setting them to NA
##################################################################################################################

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

rownames(GE) <- GE$ID
GE_SUB <- GE[,-c(1:9)]
GE_t <- as.data.frame(t(GE_SUB))

#This creates columns to merge by
GE_t$Name <- rownames(GE_t)
SUB_clean$Name <- rownames(SUB_clean)
ALL <- merge(SUB_clean, GE_t, by = "Name")

#This converts the one SNP to a numeric
ALL$S1_277994336 <- as.numeric(ALL$S1_277994336)
ALL$S4_2010582 <- as.numeric(ALL$S4_2010582)
ALL$S4_5674254 <- as.numeric(ALL$S4_5674254)
ALL$S5_190520703 <- as.numeric(ALL$S5_190520703 )
ALL$S6_175749381 <- as.numeric(ALL$S6_175749381)
ALL$S8_163979957 <- as.numeric(ALL$S8_163979957)
ALL$S8_164038154 <- as.numeric(ALL$S8_164038154)

#This drops the Name column used to merge
ALL <- ALL[,-c(1)]

#This drops some of the variables that seem non significant 
model5 <- glm(KMeans_K10L2 ~ + S1_277994336 + S4_2010582 + S4_5674254 + S5_190520703 + S6_175749381 + S8_163979957 + S8_164038154 + wind + vapr +  sq3_simp + sq4_simp + elev + bio_10 + bio_18 + bio_4 + bio_15 + sq1_simp + sq5_simp + sq6_simp + KMeans_BChrom + srad + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial", data = ALL)
summary(model5)

#This calculates the percent of deviance explained
(model5$null.deviance - model5$deviance) / model5$null.deviance * 100
#This accounts for 14.93% of deviance

#I am testing the assumptions
model5_simulationOutput <- simulateResiduals(fittedModel = model5, plot = F)

png("K10L2_AllSNPs_DHARMa_Plots_Start.png")
plot(model5_simulationOutput)
dev.off()
#There are no issues here. The residual deviance is lower than the residual degrees of freedom

step(model5)
#This is robust to rearangment 

#I started with the model recommended by step and simplified further
model6 <- glm(KMeans_K10L2 ~ PC2 + PC3 + PC6 + PC9 + PC10 + vapr + sq3_simp + sq4_simp + sq5_simp + S4_2010582 + S4_5674254 + S5_190520703 + S8_164038154, family = "binomial", data = ALL)
summary(model6)
#This model is robust to rearrangement 

ALL$sq4_simp <- gsub("1", "Good", ALL$sq4_simp)
ALL$sq4_simp <- gsub("2", "Good", ALL$sq4_simp)
ALL$sq4_simp <- gsub("3", "Bad", ALL$sq4_simp)
ALL$sq4_simp <- gsub("4", "Very Bad", ALL$sq4_simp)
ALL$sq4_simp <- factor(ALL$sq4_simp, levels = c("Good", "Bad", "Very Bad"))
summary(as.factor(ALL$sq4_simp))


ALL$sq5_simp <- gsub("1", "Good", ALL$sq5_simp)
ALL$sq5_simp <- gsub("2", "Good", ALL$sq5_simp)
ALL$sq5_simp <- gsub("3", "Good", ALL$sq5_simp)
ALL$sq5_simp <- gsub("4", "Very Bad", ALL$sq5_simp)
ALL$sq5_simp <- factor(ALL$sq5_simp, levels = c("Good", "Very Bad"))
summary(as.factor(ALL$sq5_simp))

ALL$sq3_simp <- gsub("1", "Good", ALL$sq3_simp)
ALL$sq3_simp <- gsub("2", "OK", ALL$sq3_simp)
ALL$sq3_simp <- gsub("3", "Bad", ALL$sq3_simp)
ALL$sq3_simp <- gsub("4", "Bad", ALL$sq3_simp)
ALL$sq3_simp <- factor(ALL$sq3_simp, levels = c("Good", "OK", "Bad"))
summary(as.factor(ALL$sq3_simp))

#ALL$sq3_simp <- as.factor(ALL$sq3_simp)

model6 <- glm(KMeans_K10L2 ~ S4_2010582 + S4_5674254 + S5_190520703 + S8_164038154 + vapr + sq3_simp + sq4_simp + sq5_simp + PC2 + PC3 + PC6 + PC9 + PC10, family = "binomial", data = ALL)
summary(model6)

#This calculates the percent of deviance explained
(model6$null.deviance - model6$deviance) / model6$null.deviance * 100
#This accounts for 14.36% of deviance 

#I am testing the assumptions
model6_simulationOutput <- simulateResiduals(fittedModel = model6, plot = F)

png("K10L2_AllSNPs_DHARMa_Plots_End.png")
plot(model6_simulationOutput)
dev.off()
#There are no issues here. The residual deviance is lower than the residual degrees of freedom


#This writes out the All data frame in case I need it later
write.csv(ALL, "~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/K10L2_All_Variables.csv", row.names = FALSE, quote = FALSE)

ALL <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/K10L2_All_Variables.csv")

###################################################################################
#This saves the model coefficients for plotting
###################################################################################

#This extracts the variable names
Name <- names(coef(model6))[-c(1)]
#This extracts the beta values
Beta <- coef(model6)[-c(1)]
#This extracts the p values
p <- coef(summary(model6))[,4][-c(1)]

SUM <- data.frame("Variable"=Name, "Beta"=Beta, "p"=p)
SUM$CDH <- "K10L2"

write.csv(SUM, "K10L2_GLM_FinalModelSNP.csv")

###################################################################################
#This does variance partitioning
###################################################################################

#In order for this to work I need to drop anything with missing data. I loose over half the data this way
ALL_Comp <- na.omit(ALL) 

#Several Variables drop out I removed them stepwise based on significance 
modelFULL <- glm(KMeans_K10L2 ~ S4_2010582 + S4_5674254 + S5_190520703 + S8_164038154 + vapr + sq3_simp + sq4_simp + sq5_simp + PC2 + PC3 + PC6 + PC9 + PC10, family = "binomial", data = ALL_Comp)
summary(modelFULL)
#This calculates the percent of deviance explained
(modelFULL$null.deviance - modelFULL$deviance) / modelFULL$null.deviance * 100


modelPC <- glm(KMeans_K10L2 ~ PC2 + PC3 + PC6 + PC9 + PC10, family = "binomial", data = ALL_Comp)
summary(modelPC)
#This calculates the percent of deviance explained
(modelPC$null.deviance - modelPC$deviance) / modelPC$null.deviance * 100


modelEV <- glm(KMeans_K10L2 ~ vapr + sq3_simp + sq4_simp + sq5_simp, family = "binomial", data = ALL_Comp)
summary(modelEV)
#This calculates the percent of deviance explained
(modelEV$null.deviance - modelEV$deviance) / modelEV$null.deviance * 100


modelSNP <- glm(KMeans_K10L2 ~ S4_2010582 + S4_5674254 + S5_190520703 + S8_164038154, family = "binomial", data = ALL_Comp)
summary(modelSNP)
#This calculates the percent of deviance explained
(modelSNP$null.deviance - modelSNP$deviance) / modelSNP$null.deviance * 100

#This creates a report on the variance partitioning
Report <- data.frame("CDH" = "K10L2","Full_Model"= (modelFULL$null.deviance - modelFULL$deviance) / modelFULL$null.deviance * 100, "Pop_Structure"= (modelPC$null.deviance - modelPC$deviance) / modelPC$null.deviance * 100, "Env"= (modelEV$null.deviance - modelEV$deviance) / modelEV$null.deviance * 100, "SNPs"=(modelSNP$null.deviance - modelSNP$deviance) / modelSNP$null.deviance * 100)

#This reads in the report
Report_Full <- read.csv("CDH_DeviancePartitioning.csv")

#This binds the Ab10 sample with rest of the report
Report_Full <- rbind(Report_Full, Report)

#This writes out the full report
write.csv(Report_Full, "CDH_DeviancePartitioning.csv", row.names = FALSE, quote = FALSE)



