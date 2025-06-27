library(ggplot2)
library(reshape2)
library(viridis)

VAR <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/CDH_DeviancePartitioning.csv")


VAR[4,1] <- "B Chr Pseudo\nCopy Number"

#This melts the data for plotting
MELT <- melt(VAR)

#This substitutes na for 0
MELT$value <- ifelse(is.na(MELT$value), 0, MELT$value)

labels <- c("Full_Model" = "Full Model", "Pop_Structure" = "Population Structure",  "Env"="Environment", "SNPs"="Genetic Modifiers" )

MELT$CDH <- factor(MELT$CDH, levels=c("Ab10", "B Chr", "K10L2", "B Chr Pseudo\nCopy Number"))

pdf("Deviance_Partitioning_Summary.pdf", height=8, width=4)
a <- ggplot(MELT, aes(x=variable, y=value, fill=variable)) +
  geom_col(position = "dodge" ) +
  scale_fill_viridis_d(option="turbo") +
  facet_wrap(~CDH, ncol=2, nrow=2) +
  scale_x_discrete(labels=labels)+
  labs(x="Variable Partition", y="Percent Deviance Explained") +
  theme_classic() +
  theme(axis.text = element_text(size=14), axis.title.y = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14), title = element_text(size=14),  axis.text.x = element_text(angle=90), axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(hjust=0.5), strip.text = element_text(size=14))
a
dev.off()

