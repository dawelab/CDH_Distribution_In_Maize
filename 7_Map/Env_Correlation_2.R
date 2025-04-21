library("ggplot2")
library("ggpubr")
library(corrplot)

#This loads in the 540 Genotyped samples as well as their Ab10 status as determined yb GBS and the WorldClim2 data
ENV <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

colnames(ENV) <- gsub("wc2.1_30s_", "", colnames(ENV))
#This subsets out data without GPS coordinates (ie Inbreds)
ENV_filt <- subset(ENV, ENV$Latitude != "" & is.na(ENV$bio_9) != "TRUE")

#This subsets to all biological variable, all soil variables, and growing season solar radiation, wind, and water vaport pressure variables
VAR <- c("srad_03", "srad_04", "srad_05", "srad_06", "srad_07", "srad_08", "srad_09", "srad_10", "wind_03", "wind_04", "wind_05", "wind_06", "wind_07", "wind_08", "wind_09", "wind_10", "vapr_03", "vapr_04", "vapr_05", "vapr_06", "vapr_07", "vapr_08", "vapr_09", "vapr_10", "bio_1", "bio_10", "bio_11", "bio_12", "bio_13", "bio_14", "bio_15", "bio_16", "bio_17", "bio_18", "bio_19", "bio_2", "bio_3", "bio_4", "bio_5", "bio_6", "bio_7", "bio_8", "bio_9", "elev", "sq1", "sq2", "sq3", "sq4", "sq5", "sq6", "sq7")
ENV_sub <- ENV_filt[,c(VAR)]

#This makes the correlation martix
CORMAT <- cor(ENV_sub)
#This makes a decent plot of the correlation matrix
pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Environmental_Variable_Correlation_Plot.pdf")
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
  SUB <- c("srad_03", "srad_07", "srad_10", "wind_06", "bio_10", "bio_12", "bio_15", "bio_18", "bio_19", "bio_2", "bio_4", "sq1", "sq3", "sq4", "sq5", "sq6")

ENV_sub <- ENV_filt[,c(SUB)]
#This makes the correlation martix
CORMAT_SUB <- cor(ENV_sub)
#This sets the lower triangle of the matrix to NA to remove redundant information
CORMAT_SUB[lower.tri(CORMAT_SUB)] <- NA
#This makes a decent plot of the correlation matrix
pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Environmental_Variable_Correlation_Plot_Sub.pdf")
corrplot(CORMAT_SUB)
dev.off()
#This melts the matrix for easy filtering
MELT_CORMAT_SUB <- melt(CORMAT_SUB, varnames = names(dimnames(CORMAT_SUB)), na.rm = FALSE, as.is = FALSE, value.name = "value")
#We will filter to remove one of every paur or variables that had >0.7 correlation
MELT_CORMAT_SUB_FILT <- subset(MELT_CORMAT_SUB, (value >= 0.7 | value <= -0.7) & value != 1)

