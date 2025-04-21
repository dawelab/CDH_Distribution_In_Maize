library(ggplot2)
library(tidyverse)


#This loades in Ab10 and Ab10 type calls 
GROUPS <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2Complete.csv")
BChrom <- vroom::vroom("BChromosome_Calls_final.csv")
# I determined that any sample called the same thing 95% of the time was

#This assigns class
BChrom_Pos <- subset(BChrom, Prop >= 0.95)
BChrom_Pos$KMeans_BChrom = "Yes"

BChrom_Neg <- subset(BChrom, Prop <= 0.05)
BChrom_Neg$KMeans_BChrom = "No"

BChrom_Ambig <- subset(BChrom, Prop > 0.05 & Prop < 0.95)
BChrom_Ambig$KMeans_BChrom = "Ambiguous"

BChrom <- rbind(BChrom_Pos, BChrom_Neg, BChrom_Ambig)

BChrom_small <- BChrom[,c("Name", "KMeans_BChrom", "Model")]
colnames(KMeans_BChrom) <- c("Name", "KMeans_BChrom", "BChrom_Model")

#This merges BChrom
GROUPS_2 <- merge(GROUPS, BChrom_small, by = "Name", all.x = TRUE)

#This writes out the final file
write.csv(GROUPS_2, file="/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromComplete.csv", row.names = FALSE, quote = FALSE)

GROUPS_3 <- subset(GROUPS_2, Data_Source != "Dawe_Lab_1" & Data_Source != "Dawe_Lab_2")
BChrom_Counts <- GROUPS_3 %>% 
  group_by(Maize_Type) %>% 
  count(KMeans_BChrom)
write.csv(BChrom_Counts, file="BChrom_Counts.csv")

