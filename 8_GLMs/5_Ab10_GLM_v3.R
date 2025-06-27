#install.packages("glmnet")
library(glmnet)
library(vcfR)
library(stringr)
library(DHARMa)

setwd("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map")

#Load the data with all of the CDH calls
GROUPS <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumCompleteWholeGenomePCAb10.csv")
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
DATA$KMeans_Ab10 <- gsub("Ambiguous", NA, DATA$KMeans_Ab10)
DATA$KMeans_Ab10 <- gsub("Ab10", 1, DATA$KMeans_Ab10)
DATA$KMeans_Ab10 <- gsub("N10", 0, DATA$KMeans_Ab10)
DATA$KMeans_Ab10 <- as.factor(DATA$KMeans_Ab10)

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

#This drops any variable with a missing Ab10 value
DATA <- subset(DATA, is.na(DATA$KMeans_Ab10) == FALSE)

#This drops everything but the Ab10 yes no and the variables
SUB <- DATA[,-c(1:14, 16:17, 19:20, 31:115, 135:141)]
#Set the rownames
rownames(SUB) <- DATA$Name
#This drops any row with missing data
SUB_clean <- na.omit(SUB)

#This tests for a relationship between Ab10 and elevation
model1 <- glm(KMeans_Ab10 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +  elev, family = "binomial", data=SUB_clean)
summary(model1)
#This calculates the percent of deiance explained
(model1$null.deviance - model1$deviance) / model1$null.deviance * 100

#I am testing the assumptions
model1_simulationOutput <- simulateResiduals(fittedModel = model1, plot = F)

png("Ab10_Elevation_DHARMa_Plots_Start.png")
plot(model1_simulationOutput)
dev.off()
#There are no issues here and the null deviance is lower than the degrees of freedom

###################################################################################

#This runs a traditional logisitic regression
model2 <- glm(KMeans_Ab10 ~ wind + vapr +  sq3_simp + sq4_simp + elev + bio_10 + bio_18 + bio_4 + bio_15 + sq1_simp + sq5_simp + sq6_simp + KMeans_BChrom + srad + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial", data = SUB_clean)

#This summarizes the model
summary(model2)

#This calculates the percent of deviance explained
(model2$null.deviance - model2$deviance) / model2$null.deviance * 100

#I am testing the assumptions
model2_simulationOutput <- simulateResiduals(fittedModel = model2, plot = F)

png("Ab10_All_DHARMa_Plots_Start.png")
plot(model2_simulationOutput)
dev.off()
#There are no issues here and the residual deviance is lower than the resiual df 


#This simplifies the model based on model AIC
step(model2)
#This model is robust to arrangment

#I extracted the model provided by step and then simplified further to get the below
model3 <- glm(KMeans_Ab10 ~ PC1 + PC8 + PC9 + sq1_simp, family = "binomial", data = SUB_clean)
summary(model3)
#This is robust to rearrangment during simplification


#This simplifies the model further by collapsing the soil variables
SUB_clean$sq1_simp <- gsub("1", "Good", SUB_clean$sq1_simp)
SUB_clean$sq1_simp <- gsub("2", "Good", SUB_clean$sq1_simp)
SUB_clean$sq1_simp <- gsub("3", "Good", SUB_clean$sq1_simp)
SUB_clean$sq1_simp <- gsub("4", "Bad", SUB_clean$sq1_simp)
SUB_clean$sq1_simp <- factor(SUB_clean$sq1_simp, levels = c("Good", "Bad"))
summary(as.factor(SUB_clean$sq1_simp))

#This drops some of the variables that seem non significant 
model4 <- glm(KMeans_Ab10 ~ PC1 + PC8 + PC9 + sq1_simp, family = "binomial", data = SUB_clean)
summary(model4)

#This calculates the percent of deiance explained
(model4$null.deviance - model4$deviance) / model4$null.deviance * 100
#This accounts for 3.37% of deviance 

#I am testing the assumptions
model4_simulationOutput <- simulateResiduals(fittedModel = model4, plot = F)

png("Ab10_All_DHARMa_Plots_End.png")
plot(model4_simulationOutput)
dev.off()
#There are no issues here and the residual deviance is lower than the resiual df 


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
SUM$CDH <- "Ab10"

write.csv(SUM, "Ab10_GLM_FinalModel.csv")

###################################################################################
#This attempts to add the associated loci
###################################################################################

#This loads the association data from the GWAS
Ab10_assoc <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10_results_PCA.assoc.logistic.filt")

#This filters the variants to only those that are HIGHLY significant
SIG <- subset(Ab10_assoc, neglog >= 7.3)

#This writes out the significant SNPs
write.table(SIG, "Ab10_SigAssocSNP.txt", quote = FALSE, row.names = FALSE)

VCF <- vcfR::read.vcfR("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/GWAS_Mo17_Ab10.chr.landrace.filt.imputed.vcf.gz")

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

#This reclassify the soild variables as factors
ALL$sq1_simp <- as.factor(ALL$sq1_simp)
ALL$sq2_simp <- as.factor(ALL$sq2_simp)
ALL$sq3_simp <- as.factor(ALL$sq3_simp)
ALL$sq4_simp <- as.factor(ALL$sq4_simp)
ALL$sq5_simp <- as.factor(ALL$sq5_simp)
ALL$sq6_simp <- as.factor(ALL$sq6_simp)

#this sets 0 as a factor level for SNPs
ALL$S10_71198303 <- factor(ALL$S10_71198303, levels=c(0,1,2,3,4,5))
ALL$S10_71197305 <- factor(ALL$S10_71197305, levels=c(0,1,2,3,4,5))
ALL$S9_29026173 <- factor(ALL$S9_29026173, levels=c(0,1,2,3,4,5))
ALL$S9_29026157 <- factor(ALL$S9_29026157, levels=c(0,1,2,3,4,5))
ALL$S8_73817047 <- factor(ALL$S8_73817047, levels=c(0,1,2,3,4,5))
ALL$S4_193239658 <- factor(ALL$S4_193239658, levels=c(0,1,2,3,4,5))
ALL$S4_193239646 <- factor(ALL$S4_193239646, levels=c(0,1,2,3,4,5))
ALL$S4_193239638 <- factor(ALL$S4_193239638, levels=c(0,1,2,3,4,5))
ALL$S3_192562582 <- factor(ALL$S3_192562582, levels=c(0,1,2,3,4,5))
ALL$S3_192486778 <- factor(ALL$S3_192486778, levels=c(0,1,2,3,4,5))
ALL$S3_192484757 <- factor(ALL$S3_192484757, levels=c(0,1,2,3,4,5))

#This starts the model fresh 
model5 <- glm(KMeans_Ab10 ~ sq1_simp + sq5_simp + sq6_simp + KMeans_BChrom + srad + wind + vapr +  sq3_simp + sq4_simp + elev + bio_10 + bio_18 + bio_4 + bio_15 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + S3_192484757 + S3_192486778 + S3_192562582 + S4_193239638 + S4_193239646 + S4_193239658 + S8_73817047 + S10_71197305 + S10_71198303, family = "binomial", data = ALL)
summary(model5)

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


#This reshapes the data to be merged with the rest of the environmental data
rownames(GE) <- GE$ID
GE_SUB <- GE[,-c(1:9)]
GE_t <- as.data.frame(t(GE_SUB))

#This creates columns to merge by
GE_t$Name <- rownames(GE_t)
SUB_clean$Name <- rownames(SUB_clean)
ALL <- merge(SUB_clean, GE_t, by = "Name")

#This writes out the All data frame in case I need it later
write.csv(ALL, "/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10_All_Variables.csv", row.names = FALSE, quote = FALSE)

ALL <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10_All_Variables.csv")
#This drops the Name column used to merge
ALL <- ALL[,-c(1)]

#this sets codes all of the SNPs as numerics in an additive manner
ALL$S10_71198303 <- as.numeric(ALL$S10_71198303)
ALL$S10_71197305 <- as.numeric(ALL$S10_71197305)
ALL$S9_29026173 <- as.numeric(ALL$S9_29026173)
ALL$S9_29026157 <- as.numeric(ALL$S9_29026157)
ALL$S8_73817047 <- as.numeric(ALL$S8_73817047)
ALL$S4_193239658 <- as.numeric(ALL$S4_193239658)
ALL$S4_193239646 <- as.numeric(ALL$S4_193239646)
ALL$S4_193239638 <- as.numeric(ALL$S4_193239638)
ALL$S3_192562582 <- as.numeric(ALL$S3_192562582)
ALL$S3_192486778 <- as.numeric(ALL$S3_192486778)
ALL$S3_192484757 <- as.numeric(ALL$S3_192484757)

#This filters to remove all missing values
ALL <- ALL[complete.cases(ALL), ]

#This reclassify the soil variables as factors
ALL$sq1_simp <- as.factor(ALL$sq1_simp)
ALL$sq2_simp <- as.factor(ALL$sq2_simp)
ALL$sq3_simp <- as.factor(ALL$sq3_simp)
ALL$sq4_simp <- as.factor(ALL$sq4_simp)
ALL$sq5_simp <- as.factor(ALL$sq5_simp)
ALL$sq6_simp <- as.factor(ALL$sq6_simp)


#We have determined that S9_29026173 is likely an artifact from a homolog on Ab10 and are choosing to exclude it 
#This drops some of the variables that seem non significant. I removed them stepwise based on p values
model5 <- glm(KMeans_Ab10 ~ sq1_simp + sq5_simp + sq6_simp + KMeans_BChrom + srad + wind + vapr +  sq3_simp + sq4_simp + elev + bio_10 + bio_18 + bio_4 + bio_15 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + S3_192484757 + S3_192486778 + S3_192562582 + S4_193239638 + S4_193239646 + S4_193239658 + S8_73817047 + S9_29026157 + S10_71197305 + S10_71198303, family = "binomial", data = ALL)
summary(model5)

#I am testing the assumptions
model5_simulationOutput <- simulateResiduals(fittedModel = model5, plot = F)

png("Ab10_AllSNPs_DHARMa_Plots_Start.png")
plot(model5_simulationOutput)
dev.off()
#There are no issues here and the residual deviance is lower than the resiual df 

step(model5)
#This model is robust to rearrangement 

model6 <- glm(KMeans_Ab10 ~ S3_192484757 + S3_192562582 + S4_193239646 + S8_73817047 + S9_29026157 + S10_71197305 + PC4 + PC6 + PC8 + PC9, family = "binomial", data = ALL)
summary(model6)
#This model is robust to rearrangement and the residual deviance is less than the residual degrees of freedom

#This calculates the percent of deviance explained
(model6$null.deviance - model6$deviance) / model6$null.deviance * 100
#This accounts for 10.71% of deviance 

#I am testing the assumptions
model6_simulationOutput <- simulateResiduals(fittedModel = model6, plot = F)

png("Ab10_AllSNPs_DHARMa_Plots_End.png")
plot(model6_simulationOutput)
dev.off()
#There are no issues here 



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
SUM$CDH <- "Ab10"

write.csv(SUM, "Ab10_GLM_FinalModelSNP.csv")


###################################################################################
#This does variance partitioning
###################################################################################

#In order for this to work I need to drop anything with missing data. I loose over half the data this way
ALL_Comp <- na.omit(ALL) 

#Several Variables drop out I removed them stepwise based on significance 
modelFULL <- glm(KMeans_Ab10 ~ S3_192484757 + S3_192562582 + S4_193239646 + S8_73817047 + S9_29026157 + S10_71197305 + PC4 + PC6 + PC8 + PC9, family = "binomial", data = ALL_Comp)
summary(modelFULL)
#This calculates the percent of deviance explained
(modelFULL$null.deviance - modelFULL$deviance) / modelFULL$null.deviance * 100

#Several Variables drop out I removed them stepwise based on significance 
modelSNPs <- glm(KMeans_Ab10 ~ S3_192484757 + S3_192562582 + S4_193239646 + S8_73817047 + S9_29026157 + S10_71197305, family = "binomial", data = ALL_Comp)
summary(modelSNPs)
#This calculates the percent of deviance explained
(modelSNPs$null.deviance - modelSNPs$deviance) / modelSNPs$null.deviance * 100


#Several Variables drop out I removed them stepwise based on significance 
modelPOP <- glm(KMeans_Ab10 ~ PC4 + PC6 + PC8 + PC9, family = "binomial", data = ALL_Comp)
summary(modelPOP)
#This calculates the percent of deviance explained
(modelPOP$null.deviance - modelPOP$deviance) / modelPOP$null.deviance * 100

#This reads in the report
Report_Full <- read.csv("CDH_DeviancePartitioning.csv")

#This creates a report on the variance partitioning
Report <- data.frame("CDH" = "Ab10","Full_Model"= (modelFULL$null.deviance - modelFULL$deviance) / modelFULL$null.deviance * 100, "Pop_Structure"= (modelPOP$null.deviance - modelPOP$deviance) / modelPOP$null.deviance * 100
, "Env"= NA, "SNPs"= (modelSNPs$null.deviance - modelSNPs$deviance) / modelSNPs$null.deviance * 100)

#This binds the Ab10 sample with rest of the report
Report_Full <- rbind(Report_Full, Report)

#This writes out the full report
write.csv(Report_Full, "CDH_DeviancePartitioning.csv", row.names = FALSE, quote = FALSE)
