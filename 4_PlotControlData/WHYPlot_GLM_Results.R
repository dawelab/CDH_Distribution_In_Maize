Ab10 <- read.csv("Ab10_GLM_FinalModel.csv")

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map")

labels <- c("PC2", "PC3", "PC4", "PC7", "PC8", "PC9", "PC10", "Elevation", "Mean Diurnal Range", "Temp. Seasonality", "Mean Temp. Wettest 3 Months", "Precipitation Driest 3 Months", "Precipitation Wettest Month", "Precipitation Coldest 3 Months", "soil nutrient availability\nvery severly limiting", "soil nutrient availability\nvery severly limiting")

Ab10$Direction <- ifelse(ALL$Beta < 0, "Negative", "Positive")

pdf("Ab10_GLM_Graph.pdf", height=3, width=6)
a <- ggplot(Ab10, aes(x=Variable, y=neglog, size=abs(Beta), color = Direction, shape=Direction))+
     scale_shape_manual(values=c("Negative"=6, "Positive"=2)) +
     scale_color_manual(values=c("Negative"="blue", "Positive"="red")) +
     scale_x_discrete(labels=Ab10_labels) +
     geom_point() +
     scale_y_continuous(limit=c(0,10)) +
     labs(size="Absolute Value\n Beta", y="-log(p)") +
     theme(axis.text.x = element_text(angle=90), axis.title.x = element_blank(), legend.position = "right")
a
dev.off()


K10L2 <- read.csv("K10L2_GLM_FinalModel.csv")

K10L$Direction <- ifelse(K10L2$Beta < 0, "Negative", "Positive")

ALL$Variable <- factor(ALL$Variable, levels = c("PC2", "PC3", "PC4", "PC7", "PC8", "PC9", "PC10", "bio_17", "bio_8", "elev", "sq1_simp2Very Bad", "sq4_simp2Bad", "sq4_simp2Very Bad", "bio_14", "bio_13", "bio_4", "bio_19"))


BChr <- read.csv("BChr_GLM_FinalModel.csv")
CopyB <- read.csv("BChrCopyNumber_GLM_FinalModel.csv")
