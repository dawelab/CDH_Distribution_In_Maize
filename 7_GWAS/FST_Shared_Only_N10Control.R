library(StAMPP)
library(ggplot2)
library(reshape2)

#
VCF <- read.table("/Volumes/Transcend/AllData_v5.Ab10Shared.ControlsST.filt.AB.Ab10Spec.vcf", header = TRUE)

VCF_N10 <- subset(VCF, Pop == "N10")

VCF_N10$Pop <- c(rep("B542C", 3), rep("FFMM", 3), rep("W23", 3), rep("B73", 3))

VCF_FREQ <- stamppConvert(VCF_N10, type = "r")

WINDOWS <- data.frame("start" = NA, "end" = NA, "FST_FFMM-B542C"= NA, "FST_W23-B542C" =NA, "FST_B73-B542C" = NA, "FST_W23-FFMM" = NA, "FST_B73-FFMM" = NA, "P_FFMM-B542C" = NA, "P_W23-B542C" = NA, "P_B73-B542C" = NA, "P_W23-FFMM" = NA, "P_B73-FFMM" = NA, "P_B73-W23" = NA)
WINDOWS <- WINDOWS[-c(1),]  

NAME_CHOP <- colnames(VCF_FREQ [,-c(1:5)])

WINSIZE <- 20

while(length(NAME_CHOP != 0)) {
  if(length(NAME_CHOP) >= WINSIZE) {
    SUB <- cbind(VCF_FREQ[,1:5], VCF_FREQ[,NAME_CHOP[1:WINSIZE]])
  }
  if(length(NAME_CHOP) < WINSIZE) {
    SUB <- cbind(VCF_FREQ[,1:5], VCF_FREQ[,NAME_CHOP])
  }
  FST <- stamppFst(SUB, nboots=100)
  F <- as.data.frame(FST$Fsts)
  P <- as.data.frame(FST$Pvalues)
  START <- gsub("S10_", "", NAME_CHOP[1])
  END <- gsub("S10_", "", NAME_CHOP[WINSIZE])
  `FST_FFMM-B542C` <- F[2, 1]
  `FST_W23-B542C` <- F[3, 1]
  `FST_B73-B542C` <- F[4,1]
  `FST_W23-FFMM` <- F[3,2]
  `FST_B73-FFMM` <- F[4,2]
  `FST_B73-W23` <- F[4,3]
  `P_FFMM-B542C` <- P[2, 1]
  `P_W23-B542C` <- P[3, 1]
  `P_B73-B542C` <- P[4,1]
  `P_W23-FFMM` <- P[3,2]
  `P_B73-FFMM` <- P[4,2]
  `P_B73-W23` <- P[4,3]
  NEW <<- data.frame(data.frame("start" = START, "end" = END, "FST_FFMM-B542C"= `FST_FFMM-B542C`, "FST_W23-B542C" =`FST_W23-B542C`, "FST_B73-B542C" = `FST_B73-B542C`, "FST_W23-FFMM" = `FST_W23-FFMM`, "FST_B73-FFMM" = `FST_B73-FFMM`, "P_FFMM-B542C" = `P_FFMM-B542C`, "P_W23-B542C" = `P_W23-B542C`, "P_B73-B542C" = `P_B73-B542C`, "P_W23-FFMM" = `P_W23-FFMM`, "P_B73-FFMM" = `P_B73-FFMM`, "P_B73-W23" = `P_B73-W23`))
  WINDOWS <<- rbind(WINDOWS, NEW)
  SUB_NAMES <- colnames(SUB[-c(1:5)])
  if(length(NAME_CHOP) >= WINSIZE) {
    NAME_CHOP <<- NAME_CHOP[-c(1:WINSIZE/2)]
  }
  if(length(NAME_CHOP) < WINSIZE) {
    NAME_CHOP <<- NAME_CHOP[-c(1:length(NAME_CHOP))]
  }
  
}

WINDOWS$start <- as.character(WINDOWS$start)

WINDOWS_FILT <- data.frame(start = WINDOWS$start , end= WINDOWS$end)

#This filters by P value

COMP <- c("B73.B542C", "B73.FFMM" , "FFMM.B542C", "W23.B542C", "W23.FFMM")

for(i in 1:length(COMP)) {
  x <- COMP[i]
  WINDOWS_SUB <- WINDOWS[,c("start", "end", paste("FST", x, sep="_"), paste("P", x, sep="_"))]
  WINDOWS_SUB$FST_K10L2.III <- ifelse(WINDOWS_SUB[,4] >= 0.05, NA, WINDOWS_SUB[,3])
  WINDOWS_FILT <<- cbind(WINDOWS_FILT, WINDOWS_SUB[,c(3), drop = FALSE])
}



WINDOWS_MELT <- melt(WINDOWS_FILT)

WINDOWS_MELT$start <- as.numeric(WINDOWS$start)
WINDOWS_MELT$end <- as.numeric(WINDOWS$end)

#This extracts all the SNP positions
SNPs <- colnames(VCF_FREQ [,-c(1:5)])
SNPs <- gsub("S10_", "", SNPs)
SNPs <- as.data.frame(SNPs)
colnames(SNPs) <- c("Position")
SNPs$Position <- as.numeric(SNPs$Position)

pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/8_SNPs/SNP_Hist_Shared_Only_N10Control.pdf", height=1, width = 3)
ggplot(SNPs, aes(x=Position)) + 
  geom_histogram(binwidth = 1000000, color = NA, fill = "black") +
  theme_classic()+
  scale_x_continuous(breaks = c(140000000, 150000000, 160000000, 170000000, 180000000, 190000000), labels = c(140, 150, 160, 170, 180, 190), limits = c(140000000, 190000000))
dev.off()


pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/8_SNPs/FST_SharedOnly_N10Control.pdf", height = 12, width=12)
p <- ggplot(WINDOWS_MELT) +
  geom_segment(aes(x=start, xend=end, y=value, yend=value)) +
  #facet_wrap(~variable, nrow =5 , ncol = 2, labeller = as_labeller(FACET_NAMES)) +
  facet_wrap(~variable)  +
  
  #scale_x_continuous(breaks = c(140000000, 150000000, 160000000, 170000000, 180000000, 190000000), labels = c(140, 150, 160, 170, 180, 190)) +
  #scale_y_continuous(breaks=c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), labels = c("", -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)) +
  labs(x="Ab10-I (MB)", y= "FST value") +
  annotate("rect", xmin=c(140000000), xmax=c(190000000), ymin=c(0.314) , ymax=c(0.434), alpha=0.3, color=NA, fill="red") +
  annotate("rect", xmin=c(141115174,152050000, 158250000), xmax=c(142472000,156880000, 167721000), ymin=c(-0.2, -0.2, -0.2) , ymax=c(-0.1, -0.1, -0.1), alpha=1, color=NA, fill="darkgoldenrod") +
  annotate("rect", xmin=c(142472000,150656000, 157485200), xmax=c(146699300,153145000, 159356550), ymin=c(-0.2, -0.2, -0.2) , ymax=c(-0.1, -0.1, -0.1), alpha=1, color=NA, fill="deepskyblue") +
annotate("rect", xmin=c(174433450), xmax=c(182846100), ymin=c(-0.2, -0.2, -0.2) , ymax=c(-0.1, -0.1, -0.1), alpha=1, color=NA, fill="darkorange3")
p
dev.off()

