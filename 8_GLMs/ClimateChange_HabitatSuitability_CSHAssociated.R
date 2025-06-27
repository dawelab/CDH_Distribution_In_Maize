library(ggplot2)
library(reshape2)
library(viridis)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map")

TABLE <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Stand7_Area_Table.csv")

TABLE <- TABLE[,-c(2,5)]

MELT <- melt(TABLE)

FULL <- subset(MELT, Model == "Current_Full")

pdf("CurrentFullModel_CDHAssociated_SuitableArea.pdf", height=4, width=6)
g <- ggplot(FULL, aes(x=variable, y=value, color=variable)) +
  geom_point() +
  labs(y="square meters", x="CDH") +
  theme(axis.text.x=element_text(angle=90), legend.position = "none")
g
dev.off()

TEMP1 <- subset(MELT, Model == "Current_Full" & variable == "BChrom")
TEMP2 <- subset(MELT, Model == "Current_SSP" & variable == "K10L2")
TEMP3 <- rbind(TEMP1, TEMP2)
TEMP3$Model <- "Current"
TEMP4 <- subset(MELT, Model != "Current_SSP" & Model != "Current_Full")
SSP <- rbind(TEMP3, TEMP4)

SSP$value <- SSP$value/1000

pdf("SSP_CDHAssociated_SuitableArea.pdf", height=4, width=6)
g <- ggplot(SSP, aes(x=Model, y=value, color=variable)) +
  geom_point() +
  facet_wrap(~variable ) +
  scale_color_viridis_d(option="magma", begin=0.2, end=0.8) +
  labs(y="square km", x="climate change scenario") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90), legend.position = "none")
g
dev.off()

