library(ggplot2)
library(reshape2)
library(viridis)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map")

TABLE <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Stand7_Area_Table.csv")

#This drops empty rows
TABLE <- TABLE[,-c(2,5)]

#This melts the table for plotting
MELT <- melt(TABLE)

MELT$km <- MELT$value/1000

TEMP1 <- subset(MELT, Model == "Current_Full" & variable == "BChrom")
TEMP2 <- subset(MELT, Model == "Current_SSP" & variable == "K10L2")
TEMP3 <- rbind(TEMP1, TEMP2)
TEMP3$Model <- "Current"
TEMP4 <- subset(MELT, Model != "Current_SSP" & Model != "Current_Full")
MELT <- rbind(TEMP3, TEMP4)


MELT$Model <- gsub("Current", "Present\nDay", MELT$Model)
MELT$Model <- gsub("SSP126", "2°", MELT$Model)
MELT$Model <- gsub("SSP245", "3°", MELT$Model)
MELT$Model <- gsub("SSP370", "4°", MELT$Model)
MELT$Model <- gsub("SSP585", "5°", MELT$Model)
MELT$Model <- factor(MELT$Model, levels=c("Present\nDay","2°","3°","4°", "5°"))

MELT$variable <- gsub("BChrom", "B Chr.", MELT$variable)
MELT$variable <- factor(MELT$variable, levels=c("K10L2", "B Chr."))

pdf("SuitableHabitatArea.pdf", height=2, width=4)
g <- ggplot(MELT, aes(x=Model, y=km, color=variable)) +
  geom_point(size=2) +
  scale_color_viridis(discrete=TRUE, option="magma", begin=0.5, end=0.8) +
  facet_wrap(~variable)+
  labs(y="square km", x="climate change scenario") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90), legend.position = "none")
g
dev.off()

