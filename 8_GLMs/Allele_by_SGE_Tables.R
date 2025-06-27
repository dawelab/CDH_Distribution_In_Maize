library(tidyverse)

ALL <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10_All_Variables.csv")

COPY <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/6_ClassifyAllData/Ab10_ContExperimental_CopyNumber.csv")

SUB_ALL <- ALL[,c(1:2,45)]
SUB_COPY <- COPY[,c(1,2)]

MERGE <- merge(SUB_ALL, SUB_COPY, by="Name", all = TRUE)

MERGE$CopyNum <- ifelse(MERGE$Mean >= 0.04, "2", "1")

MERGE$KMeans_Ab10 <- as.factor(MERGE$KMeans_Ab10)

MERGE$S9_29026173 <- ifelse(MERGE$S9_29026173 >= 1, "Allele 2", "Allele 1")


MERGE$KMeans_Ab10 <- gsub("0", "Absent", MERGE$KMeans_Ab10 )
MERGE$KMeans_Ab10 <- gsub("1", "Present", MERGE$KMeans_Ab10 )

TABLE <- MERGE %>%
  group_by(KMeans_Ab10) %>%
  count(MERGE$S9_29026173)


TABLE <- as.data.frame(t(TABLE))

TABLE <- MERGE %>%
  group_by(KMeans_Ab10, CopyNum) %>%
  count(MERGE$S9_29026173)