library(ggplot2)
library(ggpubr)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map")

#######################################################################
#This section plots the Ab10 results without SNPs
#######################################################################

Ab10 <- read.csv("Ab10_GLM_FinalModel.csv")
Ab10 <- Ab10[,-c(1)]
#This drops anything nonsignificant. This is included for the factor variables
Ab10 <- subset(Ab10, p <= 0.01)

Ab101labels <- c("PC1"="PC1", "PC8"="PC8", "PC9"="PC9", "sq1_simpBad"="Very Severly\n Limiting\n Soil Nutrients")

Ab10$neglog <- -log10(Ab10$p)
Ab10$Direction <- ifelse(Ab10$Beta < 0, "Negative", "Positive")
Ab10$Direction <- factor(Ab10$Direction, levels=c("Positive", "Negative"))

pdf("Ab10_GLM_Graph.pdf", height=4, width=6)
a <- ggplot(Ab10, aes(x=Variable, y=neglog, size=abs(Beta), color = Direction, fill=Direction, shape=Direction))+
        scale_shape_manual(values=c("Negative"=25, "Positive"=24)) +
        scale_color_manual(values=c("Negative"="blue", "Positive"="red")) +
        scale_fill_manual(values=c("Negative"="blue", "Positive"="red")) +
     scale_size_continuous(range= c(4,10))+
     scale_x_discrete(labels=Ab101labels) +
     geom_point() +
     ylim(0,40) +
     labs(size="Absolute Value\n Beta", y="-log(p)", title="Ab10") +
        theme_classic() +
        theme(axis.text = element_text(size=7*4), axis.title.y = element_text(size=7*4), legend.text = element_text(size=7*4), legend.title = element_text(size=7*4), title = element_text(size=7*4),  axis.text.x = element_text(angle=90), axis.title.x = element_blank(), legend.position = "right", plot.title = element_text(hjust=0.5), legend.key.size = unit(15,"mm"))
a
dev.off()

#######################################################################
#This section plots the Ab10 results with SNPs
#######################################################################

Ab10 <- read.csv("Ab10_GLM_FinalModelSNP.csv")
Ab10 <- Ab10[,-c(1)]


Ab10$neglog <- -log10(Ab10$p)
Ab10$Direction <- ifelse(Ab10$Beta < 0, "Negative", "Positive")
Ab10$Direction <- factor(Ab10$Direction, levels=c("Positive", "Negative"))

#This orders the SNPs
Ab10$Variable <- factor(Ab10$Variable, levels=c("PC4", "PC6", "PC8", "PC9", "S3_192484757", "S3_192562582", "S4_193239646", "S8_73817047", "S9_29026157", "S10_71197305"))

Ab102labels <- c("S3_192484757" = "Chr3 SNP 1", "S3_192562582" = "Chr3 SNP 3", "S4_193239646" = "Chr4 SNP 2", "S8_73817047" = "Chr8 SNP 1", "S9_29026157"="Chr9 SNP 1", "S10_71197305"="Chr10 SNP 1", "PC4"="PC4", "PC6"="PC6", "PC8"="PC8", "PC9"= "PC9")

pdf("Ab10_SNPGLM_Graph.pdf", height=3, width=6)
x <- ggplot(Ab10, aes(x=Variable, y=neglog, size=abs(Beta), color = Direction, fill= Direction, shape=Direction))+
        scale_shape_manual(values=c("Negative"=25, "Positive"=24)) +
        scale_color_manual(values=c("Negative"="blue", "Positive"="red")) +
        scale_fill_manual(values=c("Negative"="blue", "Positive"="red")) +
        scale_size_continuous(range= c(4,10))+
        scale_x_discrete(labels=Ab102labels) +
        geom_point() +
        ylim(0,30) +
        labs(size="Absolute Value\n Beta", y="-log(p)", title="Ab10") +
        theme_classic() +
        theme(axis.text = element_text(size=7*4), axis.title.y = element_text(size=7*4), legend.text = element_text(size=7*4), legend.title = element_text(size=7*4), title = element_text(size=7*4),  axis.text.x = element_text(angle=90), axis.title.x = element_blank(), legend.position = "right", plot.title = element_text(hjust=0.5), legend.key.size = unit(15,"mm"))
x
dev.off()


#######################################################################
#This section plots the K10L2 results without SNPs
#######################################################################

K10L2 <- read.csv("K10L2_GLM_FinalModel.csv")
K10L2 <- K10L2[,-c(1)]
K10L2 <- subset(K10L2 , p <= 0.01)
K10L2$neglog <- -log10(K10L2$p)
K10L2$Direction <- ifelse(K10L2$Beta < 0, "Negative", "Positive")
K10L2$Direction <- factor(K10L2$Direction, levels=c("Positive", "Negative"))


K10L2$Variable <- factor(K10L2$Variable, levels = c("PC2", "PC3", "PC10", "vapr", "bio_10", "sq4_simpBad", "sq3_simpOK"))

K10L21labels <- c("PC2"="PC2", "PC3"="PC3", "PC10"="PC10", "bio_10"="Mean Temperature \n Warmest Quarter", "bio_4"="Temperature\n Seasonality", "sq4_simpBad"="Severly Limiting \n Soil Oxygen", "sq3_simpOK"="Moderatly Limiting \n Soil Rooting Conditions",  "vapr"="Water Vapor\n Pressure")


pdf("K10L2_GLM_Graph.pdf", height=4, width=6)
b <- ggplot(K10L2, aes(x=Variable, y=neglog, size=abs(Beta), color = Direction, fill=Direction, shape=Direction))+
        scale_shape_manual(values=c("Negative"=25, "Positive"=24)) +
        scale_color_manual(values=c("Negative"="blue", "Positive"="red")) +
        scale_fill_manual(values=c("Negative"="blue", "Positive"="red")) +
        scale_x_discrete(labels=K10L21labels) +
        geom_point() +
        ylim(0,40) +
        labs(size="Absolute Value\n Beta", y="-log(p)", title="K10L2") +
        theme_classic() +
        theme(axis.text = element_text(size=7*4), axis.title.y = element_text(size=7*4), legend.text = element_text(size=7*4), legend.title = element_text(size=7*4), title = element_text(size=7*4),  axis.text.x = element_text(angle=90), axis.title.x = element_blank(), legend.position = "right", plot.title = element_text(hjust=0.5), legend.key.size = unit(15,"mm"))
b
dev.off()

#######################################################################
#This section plots the K10L2 results with SNPs
#######################################################################

K10L2 <- read.csv("K10L2_GLM_FinalModelSNP.csv")
K10L2 <- K10L2[,-c(1)]
K10L2 <- subset(K10L2 , p <= 0.01)
K10L2$neglog <- -log10(K10L2$p)
K10L2$Direction <- ifelse(K10L2$Beta < 0, "Negative", "Positive")
K10L2$Direction <- factor(K10L2$Direction, levels=c("Positive", "Negative"))


K10L2$Variable <- factor(K10L2$Variable, levels = c("PC2", "PC3", "PC6", "PC9", "PC10", "S4_2010582", "S4_5674254", "S5_190520703", "S8_164038154", "sq5_simpVery Bad", "sq3_simpOK", "sq4_simpBad", "vapr"))

K10L22labels <- c("PC2"="PC2", "PC3"="PC3", "PC6"="PC6", "PC9"="PC9", "PC10"="PC10","vapr"="Water Vapor Pressure", "sq5_simpVery Bad" = "Very Severly Limiting\n Soil Salt", "sq3_simpOK" = "Moderatly Limiting\nRooting Cond.", "sq4_simpBad" = "Severly Limiting\n Soil Oxygen", "S4_2010582"="Chr4 SNP 1", "S4_5674254"="Chr4 SNP 2", "S5_190520703" = "Chr5 SNP 1", "S8_164038154"="Chr8 SNP 2")


pdf("K10L2_SNPGLM_Graph.pdf", height=3.6, width=6)
y <- ggplot(K10L2, aes(x=Variable, y=neglog, size=abs(Beta), color = Direction, fill= Direction, shape=Direction))+
        scale_shape_manual(values=c("Negative"=25, "Positive"=24)) +
        scale_color_manual(values=c("Negative"="blue", "Positive"="red")) +
        scale_fill_manual(values=c("Negative"="blue", "Positive"="red")) +
        scale_size_continuous(range= c(4,10))+
        scale_x_discrete(labels=K10L22labels) +
        geom_point() +
        ylim(0,30) +
        labs(size="Absolute Value\n Beta", y="-log(p)", title="K10L2") +
        theme_classic() +
        theme(axis.text = element_text(size=7*4), axis.title.y = element_text(size=7*4), legend.text = element_text(size=7*4), legend.title = element_text(size=7*4), title = element_text(size=7*4),  axis.text.x = element_text(angle=90), axis.title.x = element_blank(), legend.position = "right", plot.title = element_text(hjust=0.5), legend.key.size = unit(15,"mm"))
y
dev.off()


#######################################################################
#This section plots the B Chr results without SNPs
#######################################################################


BChr <- read.csv("BChr_GLM_FinalModel.csv")
BChr <- BChr[,-c(1)]
BChr <- subset(BChr , p <= 0.01)
BChr$neglog <- -log10(BChr$p)
BChr$Direction <- ifelse(BChr$Beta < 0, "Negative", "Positive")
BChr$Direction <- factor(BChr$Direction, levels=c("Positive", "Negative"))

BChr$Variable <- factor(BChr$Variable, levels = c("PC1", "PC2", "PC5", "PC8", "PC9",  "PC10", "srad", "sq5_simp3"))

BChr1labels <- c("PC1"="PC1", "PC2"="PC2", "PC5"="PC5", "PC8"="PC8", "PC9"="PC9", "PC10"="PC10", "sq5_simp3"="Severly Limiting\n Soil Salt", "srad"= "Solar Radiation")


pdf("BChr_GLM_Graph.pdf", height=4.2, width=6)
c <- ggplot(BChr, aes(x=Variable, y=neglog, size=abs(Beta), color = Direction, fill= Direction, shape=Direction))+
        scale_shape_manual(values=c("Negative"=25, "Positive"=24)) +
        scale_color_manual(values=c("Negative"="blue", "Positive"="red")) +
        scale_fill_manual(values=c("Negative"="blue", "Positive"="red")) +
        scale_x_discrete(labels=BChr1labels) +
        geom_point() +
        ylim(0,40) +
        labs(size="Absolute Value\n Beta", y="-log(p)", title=" B Chromosome") +
        theme_classic() +
        theme(axis.text = element_text(size=7*4), axis.title.y = element_text(size=7*4), legend.text = element_text(size=7*4), legend.title = element_text(size=7*4), title = element_text(size=7*4),  axis.text.x = element_text(angle=90), axis.title.x = element_blank(), legend.position = "right", plot.title = element_text(hjust=0.5), legend.key.size = unit(15,"mm"))
c
dev.off()


#######################################################################
#This section plots the B Chr results with SNPs
#######################################################################


BChr <- read.csv("BChr_GLM_FinalModelSNP.csv")
BChr <- BChr[,-c(1)]
BChr <- subset(BChr , p <= 0.01)
BChr$neglog <- -log10(BChr$p)
BChr$Direction <- ifelse(BChr$Beta < 0, "Negative", "Positive")
BChr$Direction <- factor(BChr$Direction, levels=c("Positive", "Negative"))

BChr$Variable <- factor(BChr$Variable, levels = c("PC1", "PC2", "PC5", "PC8", "PC9", "PC10", "S3_211286197", "S4_248470517", "srad"))

BChr2labels <- c("S3_211286197"="Chr3 SNP 1", "S4_248470517"="Chr4 SNP 1", "PC1"="PC1", "PC2"="PC2", "PC5"="PC5", "PC8"="PC8", "PC9"="PC9", "PC10"="PC10", "srad"="Solar\nRadiation")


pdf("BChr_SNPGLM_Graph.pdf", height=3, width=6)
z <- ggplot(BChr,  aes(x=Variable, y=neglog, size=abs(Beta), color = Direction, fill= Direction, shape=Direction))+
        scale_shape_manual(values=c("Negative"=25, "Positive"=24)) +
        scale_color_manual(values=c("Negative"="blue", "Positive"="red")) +
        scale_fill_manual(values=c("Negative"="blue", "Positive"="red")) +
        scale_size_continuous(range= c(4,10))+
        scale_x_discrete(labels=BChr2labels) +
        geom_point() +
        ylim(0,30) +
        labs(size="Absolute Value\n Beta", y="-log(p)", title=" B Chr.") +
        theme_classic() +
        theme(axis.text = element_text(size=7*4), axis.title.y = element_text(size=7*4), legend.text = element_text(size=7*4), legend.title = element_text(size=7*4), title = element_text(size=7*4),  axis.text.x = element_text(angle=90), axis.title.x = element_blank(), legend.position = "right", plot.title = element_text(hjust=0.5), legend.key.size = unit(15,"mm"))
z
dev.off()



#######################################################################
#This section plots the B chromosome copy number
#######################################################################


CopyB <- read.csv("BChrCopyNum_GLM_FinalModel.csv")
CopyB <- CopyB[,-c(1)]
CopyB <- subset(CopyB , p <= 0.01)
CopyB$neglog <- -log10(CopyB$p)
CopyB$Direction <- ifelse(CopyB$Beta < 0, "Negative", "Positive")
CopyB$Direction <- factor(CopyB$Direction, levels=c("Positive", "Negative"))

CopyB$Variable <- factor(CopyB$Variable, levels = c("PC4", "bio_4", "bio_10", "vapr"))

CopyBlabels <- c("PC4"="  PC4  ", "bio_10"="  Mean Temperature    \nWarmest Quarter   ", "bio_4"="Temperature\nSeasonality", "vapr"="Water Vapor\n Pressure")


pdf("CopyB_GLM_Graph.pdf", height=4.4, width=6)
d <- ggplot(CopyB, aes(x=Variable, y=neglog, size=abs(Beta), color = Direction, fill= Direction, shape=Direction))+
        scale_shape_manual(values=c("Negative"=25, "Positive"=24)) +
        scale_color_manual(values=c("Negative"="blue", "Positive"="red")) +
        scale_fill_manual(values=c("Negative"="blue", "Positive"="red")) +
        scale_size_continuous(range= c(4,10))+
        scale_x_discrete(labels=CopyBlabels) +
        geom_point() +
        ylim(0,30) +
        labs(size="Absolute Value\n Beta", y="-log(p)", title=" B Chr.\n Pseudo Copy Number") +
        theme_classic() +
        theme(axis.text = element_text(size=7*4), axis.title.y = element_text(size=7*4), legend.text = element_text(size=7*4), legend.title = element_text(size=7*4), title = element_text(size=7*4),  axis.text.x = element_text(angle=90), axis.title.x = element_blank(), legend.position = "right", plot.title = element_text(hjust=0.5), legend.key.size = unit(15,"mm"))
d
dev.off()

#This plots everything together

pdf("All_GLM_Graph.pdf",  height=6*2.5, width=6.*3.7)
q <- ggarrange(x,z,y,d, common.legend = TRUE, ncol=2, nrow=2)
q
dev.off()


pdf("All_NoSNPGLM_Graph.pdf",  height=6*2.5, width=6.*3.6)
q <- ggarrange(a,c,b, common.legend = TRUE, ncol=2, nrow=2)
q
dev.off()
