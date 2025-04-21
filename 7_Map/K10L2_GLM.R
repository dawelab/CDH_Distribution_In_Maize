library(glmnet)

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
DATA$KMeans_K10L2 <- gsub("Ambiguous", NA, DATA$KMeans_K10L2)
DATA$KMeans_K10L2 <- gsub("K10L2", 1, DATA$KMeans_K10L2)
DATA$KMeans_K10L2 <- gsub("N10", 0, DATA$KMeans_K10L2)
DATA$KMeans_K10L2 <- as.factor(DATA$KMeans_K10L2)
DATA$sq1 <- factor(DATA$sq5, levels = c("1", "2", "3", "4", "5", "6", "7"))
DATA$sq2 <- factor(DATA$sq5, levels = c("1", "2", "3", "4", "5", "6", "7"))
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

#This drops any variable with a missing K10L2 value
DATA <- subset(DATA, is.na(DATA$KMeans_K10L2) == FALSE)

#This tests for a relationship between K10L2 and elevation
model1 <- glm(KMeans_K10L2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +  elev, family = "binomial", data=DATA)
summary(model1)
#This calculates the percent of deiance explained
(model1$null.deviance - model1$deviance) / model1$null.deviance * 100

#These are two methods of performing variable selection. Steps is an AIC based method and leaps uses BIC which and exhaustive algorithm to search space
#This drops everything but the K10L2 yes no and the variables
SUB <- DATA[,-c(1:16, 18:20,31:115, 135:141)]
#This drops any row with missing data
SUB_clean <- na.omit(SUB)

#This separates it into a predictor matrix and a respeons matrix
pred <- as.matrix(SUB_clean[,-c(1)])
response <- as.matrix(SUB_clean[,c(1)])

# Fit LASSO logistic regression model using glmnet
lasso_model <- glmnet(pred, response, family = "binomial", alpha = 1)

# Perform cross-validation to find optimal lambda
cv_lasso <- cv.glmnet(pred, response, family = "binomial", alpha = 1)

# Get the best lambda (lambda that minimizes cross-validation error)
best_lambda <- cv_lasso$lambda.min

# Refit the model using the best lambda
lasso_best <- glmnet(pred, response, family = "binomial", alpha = 1, lambda = best_lambda)

# Extract nonzero coefficients (selected variables)
selected_variables <- rownames(coef(lasso_best))[-1][coef(lasso_best)[-1] != 0]

# Print selected variables
print(selected_variables)

#This uses the variables selected via LASSO to run a traditional logisitic regression
model2 <- glm(KMeans_K10L2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + bio_13 + bio_14 + bio_15 + bio_17 + bio_18 + bio_19 + bio_4 + bio_8 + bio_9 + elev + sq1_simp + sq2_simp + sq3_simp + sq4_simp + sq5_simp + sq6_simp, family = "binomial", data = SUB_clean)

#This calculates the percent of deviance explained
(model2$null.deviance - model2$deviance) / model2$null.deviance * 100

#This summarizes the model
summary(model2)
anova(model2)

#This drops some of the variables that seem non significant 
model3 <-  glm(KMeans_K10L2 ~ PC2 + PC3 + PC4 + PC7 + PC8 + bio_14 + bio_15 + bio_17 + bio_8 + elev + sq1_simp + sq2_simp + sq3_simp + sq4_simp + sq5_simp + sq6_simp, family = "binomial", data = SUB_clean)
summary(model3)

#This drops more variables that now seem non significant 
model4 <-  glm(KMeans_K10L2 ~ PC2 + PC3 + PC4 + PC7 + PC8 + bio_14 + bio_17 + bio_8 + elev + sq1_simp + sq4_simp, family = "binomial", data = SUB_clean)
summary(model4)

#This collapses the remaining categorical variables further
SUB_clean$sq1_simp2 <- gsub("1", "Good", SUB_clean$sq1_simp)
SUB_clean$sq1_simp2 <- gsub("2", "Good", SUB_clean$sq1_simp2)
SUB_clean$sq1_simp2 <- gsub("3", "Good", SUB_clean$sq1_simp2)
SUB_clean$sq1_simp2 <- gsub("4", "Very Bad", SUB_clean$sq1_simp2)
SUB_clean$sq1_simp2 <- factor(SUB_clean$sq1_simp2, levels = c("Good", "Very Bad"))

SUB_clean$sq4_simp2 <- gsub("1", "Good", SUB_clean$sq4_simp)
SUB_clean$sq4_simp2 <- gsub("2", "Good", SUB_clean$sq4_simp2)
SUB_clean$sq4_simp2 <- gsub("3", "Bad", SUB_clean$sq4_simp2)
SUB_clean$sq4_simp2 <- gsub("4", "Very Bad", SUB_clean$sq4_simp2)
SUB_clean$sq4_simp2 <- factor(SUB_clean$sq4_simp2, levels = c("Good", "Bad", "Very Bad"))

#This drops more variables that now seem non significant 
model5 <-  glm(KMeans_K10L2 ~ PC2 + PC3 + PC4 + PC7 + PC8 + bio_14 + bio_17 + bio_8 + elev + sq1_simp2 + sq4_simp2, family = "binomial", data = SUB_clean)
summary(model5)
#This calculates the percent of deviance explained
(model5$null.deviance - model5$deviance) / model5$null.deviance * 100

#I rearranged variables to ensure robustness
model6 <-  glm(KMeans_K10L2 ~ PC2 + PC3 + PC4 + PC7 + PC8 + bio_8 + elev + sq1_simp2 + sq4_simp2,family = "binomial", data = SUB_clean)
summary(model6)
#This calculates the percent of deviance explained
(model6$null.deviance - model6$deviance) / model6$null.deviance * 100
#This accounts for 12.05% of the deviance



#This saves the model coefficients
#This extracts the variable names
Name <- names(coef(model6))[-c(1)]
#This extracts the beta values
Beta <- coef(model6)[-c(1)]
#This extracts the p values
p <- coef(summary(model6))[,4][-c(1)]

SUM <- data.frame("Variable"=Name, "Beta"=Beta, "p"=p)
SUM$CDH <- "K10L2"

write.csv(SUM, "K10L2_GLM_FinalModel.csv")


