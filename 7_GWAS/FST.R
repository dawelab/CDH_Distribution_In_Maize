library(StAMPP)
library(ggplot2)
library(reshape2)

#
VCF <- read.table("/Volumes/Transcend/AllData_v5.Ab10Hap.HigControls.filt.AB.vcf", header = TRUE)

VCF_FREQ <- stamppConvert(VCF, type = "r")

WINDOWS <- data.frame(start = NA, end = NA, `FST_I-II` = NA, `FST_I-III` = NA, `FST_II-III` = NA, `P_I-II` = NA, `P_I-III` = NA, `P_II-III` = NA)
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
  `FST_I-II` <- F[2, 1]
  `FST_I-III` <- F[3,1]
  `FST_II-III` <- F[3,2]
  `P_I-II` <- P[2, 1]
  `P_I-III` <- P[3,1]
  `P_II-III` <- P[3,2]
  NEW <<- data.frame("start" = START, "end" = END, "FST_I-II"=`FST_I-II`, "FST_I-III" =`FST_I-III`, "FST_II-III" = `FST_II-III`, "P_I-II" =`P_I-II`, "P_I-III"= `P_I-III`, "P_II-III"=`P_II-III`) 
  WINDOWS <<- rbind(WINDOWS, NEW)
  SUB_NAMES <- colnames(SUB[-c(1:5)])
  if(length(NAME_CHOP) >= WINSIZE) {
    NAME_CHOP <<- NAME_CHOP[-c(1:WINSIZE/2)]
  }
  if(length(NAME_CHOP) < WINSIZE) {
    NAME_CHOP <<- NAME_CHOP[-c(1:length(NAME_CHOP))]
  }
  
}


#This filters by P value
WINDOWS_FILT <- data.frame(start = WINDOWS$start , end= WINDOWS$end)


COMP <- c("I.II", "I.III", "II.III")

for(i in 1:length(COMP)) {
  x <- COMP[i]
  COLNAMES <- colnames(WINDOWS_FILT)
  WINDOWS_SUB <- WINDOWS[,c("start", "end", paste("FST", x, sep="_"), paste("P", x, sep="_"))]
  WINDOWS_SUB$NEW <- ifelse(WINDOWS_SUB[,4] >= 0.05, NA, WINDOWS_SUB[,3])
  WINDOWS_FILT <<- cbind(WINDOWS_FILT, WINDOWS_SUB[,c(5), drop = FALSE])
  colnames(WINDOWS_FILT) <- c(COLNAMES, paste("FST", x, sep="_"))
}
  
WINDOWS_MELT <- melt(WINDOWS_FILT)

#This sets if there was a significant P value
WINDOWS_MELT$SIG <- ifelse(is.na(WINDOWS_MELT[,4]) == TRUE, "No", "Yes")

#This sets a dummy value for the non significant FST
WINDOWS_MELT[is.na(WINDOWS_MELT)] = 0

#This extracts all the SNP positions
SNPs <- colnames(VCF_FREQ [,-c(1:5)])
SNPs <- gsub("S10_", "", SNPs)
SNPs <- as.data.frame(SNPs)
colnames(SNPs) <- c("Position")
SNPs$Position <- as.numeric(SNPs$Position)

pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/8_SNPs/SNP_Hist.pdf", height=1, width = 3)
ggplot(SNPs, aes(x=Position)) + 
  geom_histogram(binwidth = 1000000, color = NA, fill = "black") +
  theme_classic()+
  scale_x_continuous(breaks = c(140000000, 150000000, 160000000, 170000000, 180000000, 190000000), labels = c(140, 150, 160, 170, 180, 190))
dev.off()

#This defines facet names
FACET_NAMES <- c(
  FST_I.II = "FST Between Ab10-I and Ab10-II",
  FST_I.III = "FST Between Ab10-I and Ab10-III",
  FST_II.III = "FST Between Ab10-II and Ab10-III"
)

#This sets the start and end as numeric
WINDOWS_MELT$start <- as.numeric(WINDOWS_MELT$start)
WINDOWS_MELT$end <- as.numeric(WINDOWS_MELT$end)
WINDOWS_MELT$value <- as.numeric(WINDOWS_MELT$value)

pdf("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/8_SNPs/FST.pdf", height = 4, width=3)
p <- ggplot(WINDOWS_MELT) +
  geom_segment(aes(x=start, xend=end+0.05, y=value, yend=value, color=SIG)) +
  facet_wrap(~variable, nrow = 3, ncol = 1, labeller = as_labeller(FACET_NAMES)) +
  scale_x_continuous(breaks = c(140000000, 150000000, 160000000, 170000000, 180000000, 190000000), labels = c(140, 150, 160, 170, 180, 190)) +
  scale_y_continuous(breaks=c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4), labels = c("", -0.1, 0, 0.1, 0.2, 0.3, 0.4), limits = c(-0.2, 0.5)) +
  scale_color_manual(values= c("Yes" = "black", "No"="grey70")) +
  labs(x="Ab10-I (MB)", y= "FST value", color = "Significance") +
  annotate("rect", xmin=c(140000000), xmax=c(194631421), ymin=c(0.314) , ymax=c(0.434), alpha=0.3, color=NA, fill="salmon") +
  annotate("rect", xmin=c(141115174,152050000, 158250000), xmax=c(142472000,156880000, 167721000), ymin=c(-0.2, -0.2, -0.2) , ymax=c(-0.1, -0.1, -0.1), alpha=1, color=NA, fill="darkgoldenrod") +
  annotate("rect", xmin=c(142472000,150656000, 157485200), xmax=c(146699300,153145000, 159356550), ymin=c(-0.2, -0.2, -0.2) , ymax=c(-0.1, -0.1, -0.1), alpha=1, color=NA, fill="deepskyblue") +
annotate("rect", xmin=c(174433450), xmax=c(182846100), ymin=c(-0.2, -0.2, -0.2) , ymax=c(-0.1, -0.1, -0.1), alpha=1, color=NA, fill="darkorange3") +
  theme(legend.position="bottom")
p
dev.off()


################This section calculates Fst over the entire haplotype and then the shared and specific region separatly. 
VCF <- read.table("/Volumes/Transcend/AllData_v5.Ab10Hap.HigControls.filt.AB.vcf", header = TRUE)

VCF_FREQ <- stamppConvert(VCF, type = "r")

FST_HAP <- stamppFst(VCF_FREQ, nboots=100)

FST_HAP$Fsts
FST_HAP$Pvalues

VCF <- read.table("/Volumes/Transcend/AllData_v5.Ab10Hap.HigControls.filt.AB.Ab10Spec.vcf", header = TRUE)

VCF_FREQ <- stamppConvert(VCF, type = "r")

FST_SPEC <- stamppFst(VCF_FREQ, nboots=100)

FST_SPEC$Fsts
FST_SPEC$Pvalues

VCF <- read.table("/Volumes/Transcend/AllData_v5.Ab10Hap.HigControls.filt.AB.Ab10Shared.vcf", header = TRUE)

VCF_FREQ <- stamppConvert(VCF, type = "r")

FST_SHARE<- stamppFst(VCF_FREQ, nboots=100)

FST_SHARE$Fsts
FST_SHARE$Pvalues
