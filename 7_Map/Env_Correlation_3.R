library("ggplot2")
library("ggpubr")
library(corrplot)

#This loads in the 540 Genotyped samples as well as their Ab10 status as determined yb GBS and the WorldClim2 data
ENV <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This filters to only lines with ENV_filtironmental data
colnames(ENV) <- gsub("wc2.1_30s_", "", colnames(ENV))

#This subsets out data without GPS coordinates (ie Inbreds)
ENV_filt <- subset(ENV, Latitude != "" & is.na(bio_9) != "TRUE")

#This averages the  from solar radiation, wind, and vapor pressure
srad <- ENV_filt[,c(63:74)]
srad$sradAvg <- rowMeans(srad)
srad <- as.data.frame(srad[,13])
colnames(srad) <- c("srad")

wind <- ENV_filt[,c(75:86)]
wind$windAvg <- rowMeans(wind)
wind <- as.data.frame(wind[,13])
colnames(wind) <- c("wind")

vapr <- ENV_filt[,c(87:98)]
vapr$vaprAvg <- rowMeans(vapr)
vapr <- as.data.frame(vapr[,13])
colnames(vapr) <- c("vapr")

ENV_filt <- cbind(ENV_filt, srad, wind, vapr)

#This reformats the soil data
ENV_filt$sq1 <- factor(ENV_filt$sq1, levels = c("1", "2", "3", "4", "5", "6", "7"))
ENV_filt$sq2 <- factor(ENV_filt$sq2, levels = c("1", "2", "3", "4", "5", "6", "7"))
ENV_filt$sq3 <- factor(ENV_filt$sq3, levels = c("1", "2", "3", "4", "5", "6", "7"))
ENV_filt$sq4 <- factor(ENV_filt$sq4, levels = c("1", "2", "3", "4", "5", "6", "7"))
ENV_filt$sq5 <- factor(ENV_filt$sq5, levels = c("1", "2", "3", "4", "5", "6", "7"))
ENV_filt$sq6 <- factor(ENV_filt$sq6, levels = c("1", "2", "3", "4", "5", "6", "7"))

#ENV_filt$sq1_simp <- gsub("1", "Good", ENV_filt$sq1)
#ENV_filt$sq1_simp <- gsub("2", "Good", ENV_filt$sq1_simp)
#ENV_filt$sq1_simp <- gsub("3", "Bad", ENV_filt$sq1_simp)
#ENV_filt$sq1_simp <- gsub("4", "Bad", ENV_filt$sq1_simp)
ENV_filt$sq1_simp <- gsub("5", NA, ENV_filt$sq1)
ENV_filt$sq1_simp <- gsub("6", NA, ENV_filt$sq1_simp)
ENV_filt$sq1_simp <- gsub("7", NA, ENV_filt$sq1_simp)
ENV_filt$sq1_simp <- factor(ENV_filt$sq1_simp, levels = c("1", "2", "3", "4"))

#ENV_filt$sq2_simp <- gsub("1", "Good", ENV_filt$sq2)
#ENV_filt$sq2_simp <- gsub("2", "Good", ENV_filt$sq2_simp)
#ENV_filt$sq2_simp <- gsub("3", "Bad", ENV_filt$sq2_simp)
#ENV_filt$sq2_simp <- gsub("4", "Bad", ENV_filt$sq2_simp)
ENV_filt$sq2_simp <- gsub("5", NA, ENV_filt$sq2)
ENV_filt$sq2_simp <- gsub("6", NA, ENV_filt$sq2_simp)
ENV_filt$sq2_simp <- gsub("7", NA, ENV_filt$sq2_simp)
ENV_filt$sq2_simp <- factor(ENV_filt$sq2_simp, levels = c("1", "2", "3", "4"))

#ENV_filt$sq3_simp <- gsub("1", "Good", ENV_filt$sq3)
#ENV_filt$sq3_simp <- gsub("2", "Good", ENV_filt$sq3_simp)
#ENV_filt$sq3_simp <- gsub("3", "Bad", ENV_filt$sq3_simp)
#ENV_filt$sq3_simp <- gsub("4", "Bad", ENV_filt$sq3_simp)
ENV_filt$sq3_simp <- gsub("5", NA, ENV_filt$sq3)
ENV_filt$sq3_simp <- gsub("6", NA, ENV_filt$sq3_simp)
ENV_filt$sq3_simp <- gsub("7", NA, ENV_filt$sq3_simp)
ENV_filt$sq3_simp <- factor(ENV_filt$sq3_simp, levels = c("1", "2", "3", "4"))

#ENV_filt$sq4_simp <- gsub("1", "Good", ENV_filt$sq4)
#ENV_filt$sq4_simp <- gsub("2", "Good", ENV_filt$sq4_simp)
#ENV_filt$sq4_simp <- gsub("3", "Bad", ENV_filt$sq4_simp)
#ENV_filt$sq4_simp <- gsub("4", "Bad", ENV_filt$sq4_simp)
ENV_filt$sq4_simp <- gsub("5", NA, ENV_filt$sq4)
ENV_filt$sq4_simp <- gsub("6", NA, ENV_filt$sq4_simp)
ENV_filt$sq4_simp <- gsub("7", NA, ENV_filt$sq4_simp)
ENV_filt$sq4_simp <- factor(ENV_filt$sq4_simp, levels = c("1", "2", "3", "4"))

#ENV_filt$sq5_simp <- gsub("1", "Good", ENV_filt$sq5)
#ENV_filt$sq5_simp <- gsub("2", "Good", ENV_filt$sq5_simp)
#ENV_filt$sq5_simp <- gsub("3", "Bad", ENV_filt$sq5_simp)
#ENV_filt$sq5_simp <- gsub("4", "Bad", ENV_filt$sq5_simp)
ENV_filt$sq5_simp <- gsub("5", NA, ENV_filt$sq5)
ENV_filt$sq5_simp <- gsub("6", NA, ENV_filt$sq5_simp)
ENV_filt$sq5_simp <- gsub("7", NA, ENV_filt$sq5_simp)
ENV_filt$sq5_simp <- factor(ENV_filt$sq5_simp, levels = c("1", "2", "3", "4"))

#ENV_filt$sq6_simp <- gsub("1", "Good", ENV_filt$sq6)
#ENV_filt$sq6_simp <- gsub("2", "Good", ENV_filt$sq6_simp)
#ENV_filt$sq6_simp <- gsub("3", "Bad", ENV_filt$sq6_simp)
#ENV_filt$sq6_simp <- gsub("4", "Bad", ENV_filt$sq6_simp)
ENV_filt$sq6_simp <- gsub("5", NA, ENV_filt$sq6)
ENV_filt$sq6_simp <- gsub("6", NA, ENV_filt$sq6_simp)
ENV_filt$sq6_simp <- gsub("7", NA, ENV_filt$sq6_simp)
ENV_filt$sq6_simp <- factor(ENV_filt$sq6_simp, levels = c("1", "2", "3", "4"))




#This subsets to biological variables of interest, all soil variables, and average solar radiation, wind, and water vapor pressure variables
VAR <- c("bio_10", "bio_18", "bio_4", "bio_15", "elev", "srad" ,"wind", "vapr", "sq1_simp", "sq2_simp","sq3_simp", "sq4_simp", "sq5_simp", "sq6_simp")

ENV_filt_sub <- ENV_filt[,c(VAR)]
ENV_filt_sub <- as.data.frame(sapply(ENV_filt_sub, as.numeric ))

#Drops anything that is NA
ENV_filt_sub <- na.omit(ENV_filt_sub)

#This makes the correlation martix
CORMAT <- cor(ENV_filt_sub)

#This makes a decent plot of the correlation matrix
pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/ENV_filtironmental_Variable_Correlation_Plot.pdf")
corrplot(CORMAT)
dev.off()

#This sets the lower triangle of the matrix to NA to remove redundant information
CORMAT[lower.tri(CORMAT)] <- NA
#This melts the matrix for easy filtering
MELT_CORMAT <- melt(CORMAT, varnames = names(dimnames(CORMAT)), na.rm = FALSE, as.is = FALSE, value.name = "value")
#We will filter to remove one of every paur or variables that had >0.7 correlation
MELT_CORMAT_FILT <- subset(MELT_CORMAT, (value >= 0.7 | value <= -0.7) & value != 1)


#This begins dropping variables to remove colinearity and testing for correlation
#Dropped all vapor pressure for colinearity with bio variables
#Dropped bio1, bio11, bio5, bio8, bio6, bio9, elev for colinearity with 10, maize would grow in warmest quarter seems more informative
#Dropped bio 16 and 13 for colinearity with bio 12, neither precipitation of wettest quarter or precipitation of wettest month seems particularly relevant to corn
#Dropped bio 14 and 17 in favor of bio 15, seasonality of precipitation encompasses most of what you would get from 14 and 17
#dropped 7 as it is colinear with several things and is not particularlt relevant to corn
#dropped bio 3 for bio 4
#dropped sq2 for sq1 fertilizer is unlikely to be applied in landraces
#dropped sq7 for sq3 workability is more relevant to farmer while rooting conditions is more releant to plant
#I am choosing to keep sq5(salt) and sq6(toxin) because they are both interesting to me and are not all that correlated (.78)
SUB <- c("bio_10", "bio_18", "bio_4", "bio_15", "srad" ,"wind", "vapr", "sq1_simp", "sq3_simp", "sq4_simp", "sq5_simp", "sq6_simp")

ENV_filt_sub2 <- ENV_filt_sub[,c(SUB)]
#This makes the correlation martix
CORMAT_SUB <- cor(ENV_filt_sub2)

#This makes a decent plot of the correlation matrix
pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/ENV_filtironmental_Variable_Correlation_Plot_Sub.pdf")
corrplot(CORMAT_SUB)
dev.off()

#This sets the lower triangle of the matrix to NA to remove redundant information
CORMAT_SUB[lower.tri(CORMAT_SUB)] <- NA
#This melts the matrix for easy filtering
MELT_CORMAT_SUB <- melt(CORMAT_SUB, varnames = names(dimnames(CORMAT_SUB)), na.rm = FALSE, as.is = FALSE, value.name = "value")
#We will filter to remove one of every paur or variables that had >0.7 correlation
MELT_CORMAT_SUB_FILT <- subset(MELT_CORMAT_SUB, (value >= 0.7 | value <= -0.7) & value != 1)

