library(ggplot2)
library(tidyverse)

#This loades in Ab10 and Ab10 type calls 
#from 1.7
GROUPS <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")
#from 6.2.3
Ab10 <- vroom::vroom("Kmeans_Exp_Groups_MixedControls_All.csv")
#from 6.2.7
TYPE <- vroom::vroom("RF_Ab10TypeClassifications.csv")

#This selects Ab10 positive experimental samples
Ab10 <- Ab10[,-c(127)]
Ab10$All <- rowSums(Ab10[,-c(1)])
Ab10$Perc <- Ab10$All/125

ggplot(Ab10, aes(x=Perc)) +
  geom_histogram() +
  xlab("Proportion Models Calling Ab10")
ggsave("Proportion_KMeansModels_Ab10.png")

#This assigns calss
Ab10_Pos <- subset(Ab10, Perc >= 0.95)
Ab10_Pos$KMeansClass = "Ab10"

Ab10_Neg <- subset(Ab10, Perc <= 0.05)
Ab10_Neg$KMeansClass = "N10"

Ab10_Ambig <- subset(Ab10, Perc > 0.05 & Perc < 0.95)
Ab10_Ambig$KMeansClass = "Ambiguous"

Ab10 <- rbind(Ab10_Pos, Ab10_Neg, Ab10_Ambig)

Ab10_small <- Ab10[,c("Name", "KMeansClass")]
colnames(Ab10_small) <- c("Name", "KMeans_Ab10")

#This merges Ab10 
GROUPS_2 <- merge(GROUPS, Ab10_small, by = "Name", all.x = TRUE)

#This merges Type
colnames(TYPE) <- c("Name", "RF_Ab10Type")
GROUPS_3 <- merge(GROUPS_2, TYPE, by="Name", all.x = TRUE)

#This writes out the final file
write.csv(GROUPS_3, file="Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10Complete.csv", row.names = FALSE, quote = FALSE)

GROUPS_4 <- subset(GROUPS_3, Data_Source != "Dawe_Lab_1" & Data_Source != "Dawe_Lab_2")
Chr10_Counts <- GROUPS_4 %>% 
  group_by(Maize_Type) %>% 
  count(KMeans_Ab10)
write.csv(Chr10_Counts, file="Chr10_Counts.csv", row.names = FALSE, quote = FALSE)

GROUPS_4 <- subset(GROUPS_3, Data_Source != "Dawe_Lab_1" & Data_Source != "Dawe_Lab_2")
Ab10Type_Counts <- GROUPS_4 %>% 
  group_by(Maize_Type) %>% 
  count(RF_Ab10Type)
write.csv(Ab10Type_Counts, file="Ab10Type_Counts.csv", row.names = FALSE, quote = FALSE)

Landrace_Ab10II <- subset(GROUPS_4, Maize_Type == "Landrace" & RF_Ab10Type == "Ab10-II")
write.csv(Landrace_Ab10II, file="Landrace_Ab10II.csv", row.names = FALSE, quote = FALSE)

