library(ggplot2)
library(viridis)

Ab10 <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/6_ClassifyAllData/Ab10_ContExperimental_CopyNumber.csv")

K10L2 <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/6_ClassifyAllData/K10L2_ContExperimental_CopyNumber.csv")

BChr <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/6_ClassifyAllData/BChr_ContExperimental_CopyNumber.csv")

ALL <- rbind(Ab10, K10L2, BChr)

ALL$Copy <- gsub("Unknown", "BChr\nPos.\nControl", ALL$Copy)
ALL$Copy <- gsub("Experimental", "Pos.\nExp.", ALL$Copy)

ALL$Copy <- factor(ALL$Copy, levels = c(0, 1, 2, "BChr\nPos.\nControl", "Pos.\nExp."))


ALL$CDH <- factor(ALL$CDH, levels = c("Ab10", "K10L2", "BChr"))

#I want all columns to have space for 3 bars so I am adding some blank lines
ADD <- ALL[c(1),]
ADD <- ADD[-c(1),]
ADD[1:4,] <- NA
ADD[,2] <- 0.25
ADD[,1] <- "Dummy"
ADD[1,3] <- 1
ADD[1,22] <- "BChr"
ADD[2,3] <- 2
ADD[2,22] <- "BChr"
ADD[3,3] <- "BChr\nPos.\nControl"
ADD[3,22] <- "Ab10"
ADD[4,3] <- "BChr\nPos.\nControl"
ADD[4,22] <- "K10L2"

ALL <- rbind(ALL, ADD)

pdf("CopyNumber_Paper.pdf", height=4, width=4)
a <- ggplot(ALL, aes(x=Copy, y=Mean, color=CDH, fill = CDH)) +
  scale_color_viridis(discrete=TRUE, begin= 0.1, end = 0.8, option="magma") +
  scale_fill_viridis(discrete=TRUE, begin= 0.1, end = 0.8, option="magma") +
  geom_hline(yintercept = 0.02, color="grey", linetype='dotted', linewidth=0.5) +
  geom_hline(yintercept = 0.04, color="grey", linetype='dotted', linewidth=0.5) +
  geom_hline(yintercept = 0.06, color="grey", linetype='dotted', linewidth=0.5) +
  geom_hline(yintercept = 0.08, color="grey", linetype='dotted', linewidth=0.5) +
  geom_hline(yintercept = 0.1, color="grey", linetype='dotted', linewidth=0.5) +
  geom_hline(yintercept = 0.12, color="grey", linetype='dotted', linewidth=0.5) +
  geom_hline(yintercept = 0.14, color="grey", linetype='dotted', linewidth=0.5) +
  geom_hline(yintercept = 0.16, color="grey", linetype='dotted', linewidth=0.5) +
  geom_hline(yintercept = 0.18, color="grey", linetype='dotted', linewidth=0.5) +
  geom_hline(yintercept = 0.20, color="grey", linetype='dotted', linewidth=0.5) +
  geom_boxplot(alpha=0.5) +
  labs(x="Copy Number", y="Pseudo Copy Number") +
  theme_classic() +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=12), legend.text = element_text(size=12), legend.position = "right")

a
dev.off()

