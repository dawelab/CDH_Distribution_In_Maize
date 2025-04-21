library(ggplot2)
library(tidyverse)


#This loades in Ab10 and Ab10 type calls 
GROUPS <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10Complete.csv")
K10L2 <- vroom::vroom("/Volumes/Transcend/7.6.7_K10L2_Kmeans_Exp_Groups_MixedControls_All.csv")
# I determined that any sample called the same thing 95% of the time was

#This selects Ab10 positive experimental samples
K10L2 <- K10L2[,-c(127, 128)]
K10L2$All <- rowSums(K10L2[,-c(1)])
K10L2$Perc <- K10L2$All/125

ggplot(K10L2, aes(x=Perc)) +
  geom_histogram() +
  xlab("Proportion Models Calling K10L2")
ggsave("Proportion_KMeansModels_K10L2.png")

ggplot(K10L2, aes(x=Perc)) +
  geom_histogram() +
  ylim(0,400) +
  xlab("Proportion Models Calling K10L2")
ggsave("Proportion_KMeansModels_K10L2_Zoom.png")


#This assigns calss
K10L2_Pos <- subset(K10L2, Perc >= 0.95)
K10L2_Pos$KMeansClass = "K10L2"

K10L2_Neg <- subset(K10L2, Perc <= 0.05)
K10L2_Neg$KMeansClass = "N10"

K10L2_Ambig <- subset(K10L2, Perc > 0.05 & Perc < 0.95)
K10L2_Ambig$KMeansClass = "Ambiguous"

K10L2 <- rbind(K10L2_Pos, K10L2_Neg, K10L2_Ambig)

K10L2_small <- K10L2[,c("Name", "KMeansClass")]
colnames(K10L2_small) <- c("Name", "KMeans_K10L2")

#This merges Ab10 
GROUPS_2 <- merge(GROUPS, K10L2_small, by = "Name", all.x = TRUE)

#This writes out the final file
write.csv(GROUPS_2, file="/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2Complete.csv", row.names = FALSE, quote = FALSE)

GROUPS_3 <- subset(GROUPS_2, Data_Source != "Dawe_Lab_1" & Data_Source != "Dawe_Lab_2")
K10L2_Counts <- GROUPS_3 %>% 
  group_by(Maize_Type) %>% 
  count(KMeans_K10L2)
write.csv(K10L2_Counts, file="K10L2_Counts.csv")
