library(ggplot2)
library(vcfR)
library(ggpubr)

setwd("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/8_SNPs")

########################################################
#Looks at some of the filtering
########################################################

VAR <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/plink.lmiss", sep="")
SAMP <-  read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/plink.imiss", sep="")

hist(VAR$F_MISS)
hist(SAMP$F_MISS)

MAF <- read.table("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/MAF_check.frq", header = TRUE)
png("MAF_hist.png")
hist(MAF$MAF)
dev.off()

HWE <- read.table("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/plink.hwe", header = TRUE)
hist(HWE)

#This creates a list of BChrom positive individuals
GROUPS <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/6.5_Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumComplete.csv")

BChrom <- subset(GROUPS, KMeans_BChrom == "Yes")
DF <- data.frame("FID"=BChrom$Name, "IID"=BChrom$Name)
write.table(DF, "BChrom_Positive_Only.txt", quote = FALSE, row.names = FALSE)

Land <- subset(GROUPS, Maize_Type == "Landrace")
DF <- data.frame("FID"=Land$Name, "IID"=Land$Name)
write.table(DF, "Landrace_Only.txt", quote = FALSE, row.names = FALSE)




########################################################
#This generates Manhattan plots
########################################################

Ab10_assoc <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10_results_PCA.assoc.logistic", sep="")
Ab10_assoc$neglogp <- -log10(Ab10_assoc$P)
Ab10_assoc$Direction <- ifelse(Ab10_assoc$BETA < 0, "Negative", "Positive")
#This drops any association to a PC
Ab10_assoc <- subset(Ab10_assoc, TEST == "ADD")

#Drop NA
Ab10_assoc <- subset(Ab10_assoc, is.na(BETA) == FALSE)

pdf("Ab10_GWAS.pdf", height=2, width=6)
ggplot(Ab10_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  ggtitle("Ab10") +
  geom_hline(yintercept=7.3, color="black", linetype="dotted") +
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5))
dev.off()


#This identifies the SNP name of SNPs that overlap orthologs of Ab10 genes
Ab10_ortho <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10_results_PCA.assoc.logistic.Ab10Ortho.bed", sep="")
colnames(Ab10_ortho) <- c("SNP Chr", "SNP Start", "SNP End", "Gene Chr", "Gene Start", "Gene End")

#This generates the SNP namame for those that overlap Ab10 orthologs
Ab10_ortho <- subset(Ab10_ortho, Ab10_ortho$`Gene Chr` != ".")
Ab10_ortho_SNPs <- unique(Ab10_ortho[,c(1:2)])
Ab10_ortho_SNPs$`SNP Chr` <- gsub("chr", "", Ab10_ortho_SNPs$`SNP Chr`)
Ab10_ortho_SNPs$Name <- paste0("S", Ab10_ortho_SNPs$`SNP Chr`,"_", Ab10_ortho_SNPs$`SNP Start` )

#This drops any SNPs by that name
Ab10_assoc <- Ab10_assoc[-c(which(Ab10_assoc$SNP %in% Ab10_ortho_SNPs$Name)),]


#This identifies the SNP name of SNPs that overlap TEs
Ab10_TE <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10_results_PCA.assoc.logistic.TE.bed", sep="", header=FALSE)
colnames(Ab10_TE) <- c("SNP Chr", "SNP Start", "SNP End", "TE Chr", "TE Tool",  "TE method", "TE Start", "TE End", "Misc1", "Misc2", "Misc3", "Misc4")

#This generates the SNP namame for those that overlap Ab10 orthologs
Ab10_TE <- subset(Ab10_TE, Ab10_TE$`TE Chr` != ".")
Ab10_TE_SNPs <- unique(Ab10_TE[,c(1:2)])
Ab10_TE_SNPs$`SNP Chr` <- gsub("chr", "", Ab10_TE_SNPs$`SNP Chr`)
Ab10_TE_SNPs$Name <- paste0("S", Ab10_TE_SNPs$`SNP Chr`,"_", Ab10_TE_SNPs$`SNP Start` )

#This drops any SNPs by that name
Ab10_assoc <- Ab10_assoc[-c(which(Ab10_assoc$SNP %in% Ab10_TE_SNPs$Name)),]


#This plots a histogram
pdf("Ab10_GWAS_Filt_Hist.pdf", height=2, width=6)
w <- ggplot(Ab10_assoc, aes(x=BP)) +
  geom_histogram(binwidth=1000000) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  ggtitle(paste("Ab10\n", nrow(Ab10_assoc))) +
  labs(x="Genomic Location", y="count") +
  theme_classic() +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5), legend.position = "none", axis.text.y = element_text(size=14))
w
dev.off()

#This makes the filtered manhattan plot
pdf("Ab10_GWAS_Filt.pdf", height=2, width=6)
a <- ggplot(Ab10_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  #scale_y_continuous(breaks=c(0,10,20,30,80), limits=c(0,90)) +
  ggtitle("Ab10") +
  geom_hline(yintercept=7.3, color="black", linetype="dotted") +
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5), legend.position = "none", axis.text.y = element_text(size=14))
a
dev.off()

#This isolates the significant SNPs only
Ab10_GWAS_SIG <- subset(Ab10_assoc, neglogp > 7.3)

#I blasted the area surrounding these SNPs and have here a list of all SNPs with greater than 80% identity to Ab10 and a hit of atleast 62 bp. I checked all of the K10L2 SNPs and found none.
Ab10_DROP <- c("S2_13304128", "S2_13304149", "S2_49085196" )

#This drops any SNPs by that name
Ab10_assoc <- Ab10_assoc[-c(which(Ab10_assoc$SNP %in% Ab10_DROP)),]

#This makes the filtered manhattan plot
pdf("Ab10_GWAS_Filt2.pdf", height=2, width=6)
a <- ggplot(Ab10_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  #scale_y_continuous(breaks=c(0,10,20,30,80), limits=c(0,90)) +
  ggtitle("Ab10") +
  geom_hline(yintercept=7.3, color="black", linetype="dotted") +
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  scale_y_continuous(limits=c(0,50), breaks=c(0,10,20,30,40,50)) +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5), legend.position = "none", axis.text.y = element_text(size=14))
a
dev.off()

write.csv(Ab10_assoc, "~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10_results_PCA.assoc.logistic.filt", row.names = FALSE, quote = FALSE)

#################################################################################

K10L2_assoc <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/K10L2_results_PCA.assoc.logistic", sep="")
K10L2_assoc$neglogp <- -log10(K10L2_assoc$P)
K10L2_assoc$Direction <- ifelse(K10L2_assoc$BETA < 0, "Negative", "Positive")

pdf("K10L2_GWAS.pdf", height=2, width=6)
ggplot(K10L2_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  ggtitle("K10L2") +
  geom_hline(yintercept=7.3, color="grey80", linetype="dotted") +
  
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5))
dev.off()

#This identifies the SNP name of SNPs that overlap orthologs of K10L2 genes
K10L2_ortho <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/K10L2_results_PCA.assoc.logistic.K10L2Ortho.bed", sep="")
colnames(K10L2_ortho) <- c("SNP Chr", "SNP Start", "SNP End", "Gene Chr", "Gene Start", "Gene End")

#This generates the SNP namame for those that overlap K10L2 orthologs
K10L2_ortho <- subset(K10L2_ortho, K10L2_ortho$`Gene Chr` != ".")
K10L2_ortho_SNPs <- unique(K10L2_ortho[,c(1:2)])
K10L2_ortho_SNPs$`SNP Chr` <- gsub("chr", "", K10L2_ortho_SNPs$`SNP Chr`)
K10L2_ortho_SNPs$Name <- paste0("S", K10L2_ortho_SNPs$`SNP Chr`,"_", K10L2_ortho_SNPs$`SNP Start` )

#This drops any SNPs by that name
K10L2_assoc <- K10L2_assoc[-c(which(K10L2_assoc$SNP %in% K10L2_ortho_SNPs$Name)),]

#This identifies the SNP name of SNPs that overlap TEs
K10L2_TE <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/K10L2_results_PCA.assoc.logistic.TE.bed", sep="", header=FALSE)
colnames(K10L2_TE) <- c("SNP Chr", "SNP Start", "SNP End", "TE Chr", "TE Tool",  "TE method", "TE Start", "TE End", "Misc1", "Misc2", "Misc3", "Misc4")

#This generates the SNP namame for those that overlap K10L2 orthologs
K10L2_TE <- subset(K10L2_TE, K10L2_TE$`TE Chr` != ".")
K10L2_TE_SNPs <- unique(K10L2_TE[,c(1:2)])
K10L2_TE_SNPs$`SNP Chr` <- gsub("chr", "", K10L2_TE_SNPs$`SNP Chr`)
K10L2_TE_SNPs$Name <- paste0("S", K10L2_TE_SNPs$`SNP Chr`,"_", K10L2_TE_SNPs$`SNP Start` )

#This drops any SNPs by that name
K10L2_assoc <- K10L2_assoc[-c(which(K10L2_assoc$SNP %in% K10L2_TE_SNPs$Name)),]

#This drops any association to a PC
K10L2_assoc <- subset(K10L2_assoc, TEST == "ADD")

#Drop NA
K10L2_assoc <- subset(K10L2_assoc, is.na(BETA) == FALSE)

#This plots a histogram
pdf("K10L2_GWAS_Filt_Hist.pdf", height=2, width=6)
x <- ggplot(K10L2_assoc, aes(x=BP)) +
  geom_histogram(binwidth=1000000) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  ggtitle(paste("K10L2\n", nrow(K10L2_assoc))) +
  labs(x="Genomic Location", y="count") +
  theme_classic() +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5), legend.position = "none", axis.text.y = element_text(size=14))
x
dev.off()

pdf("K10L2_GWAS_Filt.pdf", height=2, width=6)
b <- ggplot(K10L2_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  ggtitle("K10L2") +
  geom_hline(yintercept=7.3, color="black", linetype="dotted") +
  
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  #ylim(0,30) +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5), legend.position = "none", axis.text.y = element_text(size=14))
b
dev.off()

#This isolates the significant SNPs only
K10L2_GWAS_SIG <- subset(K10L2_assoc, neglogp > 7.3)

#I blasted the area surrounding these SNPs and have here a list of all SNPs with greater than 80% identity to K10L2. I checked all of the K10L2 SNPs and found none.
#K10L2_DROP <- c("")

#This drops any SNPs by that name
#K10L2_assoc <- K10L2_assoc[-c(which(K10L2_assoc$SNP %in% K10L2_DROP)),]

pdf("K10L2_GWAS_Filt2.pdf", height=2, width=6)
b <- ggplot(K10L2_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  ggtitle("K10L2") +
  geom_hline(yintercept=7.3, color="black", linetype="dotted") +
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  scale_y_continuous(limits=c(0,50), breaks=c(0,10,20,30,40,50)) +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5), legend.position = "none", axis.text.y = element_text(size=14))
b
dev.off()

write.csv(K10L2_assoc, "~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/K10L2_results_PCA.assoc.logistic.filt", row.names = FALSE, quote = FALSE)

#################################################################################

BChrom_assoc <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChrom_results_PCA.assoc.logistic", sep="")
BChrom_assoc$neglogp <- -log10(BChrom_assoc$P)
BChrom_assoc$Direction <- ifelse(BChrom_assoc$BETA < 0, "Negative", "Positive")

pdf("BChrom_GWAS.pdf", height=2, width=6)
ggplot(BChrom_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  ggtitle("B Chromosome") +
  geom_hline(yintercept=7.3, color="grey80", linetype="dotted") +
  
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5))
dev.off()

#This identifies the SNP name of SNPs that overlap orthologs of B chromosomes genes
BChrom_ortho <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChrom_results_PCA.assoc.logistic.BChromOrtho.bed", sep="")
colnames(BChrom_ortho) <- c("SNP Chr", "SNP Start", "SNP End", "Gene Chr", "Gene Start", "Gene End")

#This generates the SNP namame for those that overlap BChrom orthologs
BChrom_ortho <- subset(BChrom_ortho, BChrom_ortho$`Gene Chr` != ".")
BChrom_ortho_SNPs <- unique(BChrom_ortho[,c(1:2)])
BChrom_ortho_SNPs$`SNP Chr` <- gsub("chr", "", BChrom_ortho_SNPs$`SNP Chr`)
BChrom_ortho_SNPs$Name <- paste0("S", BChrom_ortho_SNPs$`SNP Chr`,"_", BChrom_ortho_SNPs$`SNP Start` )

#This drops any SNPs by that name
BChrom_assoc <- BChrom_assoc[-c(which(BChrom_assoc$SNP %in% BChrom_ortho_SNPs$Name)),]

#This identifies the SNP name of SNPs that overlap TEs
BChrom_TE <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChrom_results_PCA.assoc.logistic.TE.bed", sep="", header=FALSE)
colnames(BChrom_TE) <- c("SNP Chr", "SNP Start", "SNP End", "TE Chr", "TE Tool",  "TE method", "TE Start", "TE End", "Misc1", "Misc2", "Misc3", "Misc4")

#This generates the SNP namame for those that overlap BChrom orthologs
BChrom_TE <- subset(BChrom_TE, BChrom_TE$`TE Chr` != ".")
BChrom_TE_SNPs <- unique(BChrom_TE[,c(1:2)])
BChrom_TE_SNPs$`SNP Chr` <- gsub("chr", "", BChrom_TE_SNPs$`SNP Chr`)
BChrom_TE_SNPs$Name <- paste0("S", BChrom_TE_SNPs$`SNP Chr`,"_", BChrom_TE_SNPs$`SNP Start` )

#This drops any SNPs by that name
BChrom_assoc <- BChrom_assoc[-c(which(BChrom_assoc$SNP %in% BChrom_TE_SNPs$Name)),]


#This drops any association to a PC
BChrom_assoc <- subset(BChrom_assoc, TEST == "ADD")

#Drop NA
BChrom_assoc <- subset(BChrom_assoc, is.na(BETA) == FALSE)


#This plots a histogram

pdf("BChrom_GWAS_Filt_Hist.pdf", height=2, width=6)
y <- ggplot(BChrom_assoc, aes(x=BP)) +
  geom_histogram(binwidth=1000000) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  ggtitle(paste("B Chr.\n", nrow(BChrom_assoc))) +
  labs(x="Genomic Location", y="count") +
  theme_classic() +
  #scale_y_continuous(limits = c(0,60), breaks=c(0,20,40,60)) +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5), legend.position = "none", axis.text.y = element_text(size=14))
y
dev.off()

pdf("BChrom_GWAS_Filt.pdf", height=2, width=6)
c <- ggplot(BChrom_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  ggtitle("B Chr.") +
  geom_hline(yintercept=7.3, color="black", linetype="dotted") +
  
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  #ylim(0,30) +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5), legend.position = "none", axis.text.y = element_text(size=14))
c
dev.off()

#This isolates the significant SNPs only
BChrom_GWAS_SIG <- subset(BChrom_assoc, neglogp > 7.3)

#I blasted the area surrounding these SNPs and have here a list of all SNPs with greater than 80% identity to the B 
DROP <- c("S2_90555553", 	"S2_90555560", "S3_5366715", "S3_5366738", "S6_124809191", "S6_124809195", "S8_199282680", "S6_124809146", "S4_178662947", "S6_124809142")

#This drops any SNPs by that name
BChrom_assoc <- BChrom_assoc[-c(which(BChrom_assoc$SNP %in% DROP)),]


#This one is on a scaffold S3_5366738, S3_5366715

pdf("BChrom_GWAS_Filt2.pdf", height=2, width=6)
c <- ggplot(BChrom_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  ggtitle("B Chr.") +
  geom_hline(yintercept=7.3, color="black", linetype="dotted") +
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  scale_y_continuous(limits=c(0,50), breaks=c(0,10,20,30,40,50)) +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5), legend.position = "none", axis.text.y = element_text(size=14))
c
dev.off()

write.csv(BChrom_assoc, "~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChrom_results_PCA.assoc.logistic.filt", row.names = FALSE, quote = FALSE)

#################################################################################

CopyBChrom_assoc <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChromCopyNum_results_PCA.assoc.linear", sep="")
CopyBChrom_assoc$neglogp <- -log10(CopyBChrom_assoc$P)
CopyBChrom_assoc$Direction <- ifelse(CopyBChrom_assoc$BETA < 0, "Negative", "Positive")

#Select only the correct tests
CopyBChrom_assoc <- subset(CopyBChrom_assoc, TEST == "ADD")

#Drop NA
CopyBChrom_assoc <- subset(CopyBChrom_assoc, is.na(BETA) == FALSE)


pdf("CopyBChrom_GWAS.pdf", height=2, width=6)
d <- ggplot(CopyBChrom_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  ggtitle("B Chr. Pseudo Copy Number") +
  geom_hline(yintercept=7.3, color="black", linetype="dotted") +
  
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  #ylim(0,30) +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5), legend.position = "none", axis.text.y = element_text(size=14))
d
dev.off()

#This drops any SNPs that are B chrom homologs
CopyBChrom_ortho <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChromCopyNum_results_PCA.assoc.linear.BChromOrtho.bed", sep="")
colnames(CopyBChrom_ortho) <- c("SNP Chr", "SNP Start", "SNP End", "Gene Chr", "Gene Start", "Gene End")

#This generates the SNP namame for those that overlap CopyBChrom orthologs
#There are none so this fails. I commented it out
CopyBChrom_ortho <- subset(CopyBChrom_ortho, CopyBChrom_ortho$`Gene Chr` != ".")
CopyBChrom_ortho_SNPs <- unique(CopyBChrom_ortho[,c(1:2)])
CopyBChrom_ortho_SNPs$`SNP Chr` <- gsub("chr", "", CopyBChrom_ortho_SNPs$`SNP Chr`)
CopyBChrom_ortho_SNPs$Name <- paste0("S", CopyBChrom_ortho_SNPs$`SNP Chr`,"_", CopyBChrom_ortho_SNPs$`SNP Start` )

#This drops any SNPs by that name
CopyBChrom_assoc <- CopyBChrom_assoc[-c(which(CopyBChrom_assoc$SNP %in% CopyBChrom_ortho_SNPs$Name)),]

#This identifies the SNP name of SNPs that overlap TEs
CopyBChrom_TE <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChromCopyNum_results_PCA.assoc.linear.TE.bed", sep="", header=FALSE)
colnames(CopyBChrom_TE) <- c("SNP Chr", "SNP Start", "SNP End", "TE Chr", "TE Tool",  "TE method", "TE Start", "TE End", "Misc1", "Misc2", "Misc3", "Misc4")

#This generates the SNP namame for those that overlap CopyBChrom orthologs
CopyBChrom_TE <- subset(CopyBChrom_TE, CopyBChrom_TE$`TE Chr` != ".")
CopyBChrom_TE_SNPs <- unique(CopyBChrom_TE[,c(1:2)])
CopyBChrom_TE_SNPs$`SNP Chr` <- gsub("chr", "", CopyBChrom_TE_SNPs$`SNP Chr`)
CopyBChrom_TE_SNPs$Name <- paste0("S", CopyBChrom_TE_SNPs$`SNP Chr`,"_", CopyBChrom_TE_SNPs$`SNP Start` )

#This drops any SNPs by that name
CopyBChrom_assoc <- CopyBChrom_assoc[-c(which(CopyBChrom_assoc$SNP %in% CopyBChrom_TE_SNPs$Name)),]


pdf("CopyBChrom_GWAS_Hist.pdf", height=2, width=6)
z <- ggplot(CopyBChrom_assoc, aes(x=BP)) +
  geom_histogram(binwidth=1000000) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  ggtitle(paste("B Chr. Pseudo Copy Number\n", nrow(CopyBChrom_assoc))) +
  labs(x="Genomic Location", y="count") +
  theme_classic() +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5), legend.position = "none", axis.text.y = element_text(size=14))
z
dev.off()

pdf("CopyBChrom_GWAS_Filt.pdf", height=2, width=6)
d <- ggplot(CopyBChrom_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  ggtitle("B Chr. Pseudo Copy Number") +
  geom_hline(yintercept=7.3, color="black", linetype="dotted") +
  
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  #ylim(0,30) +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5), legend.position = "none", axis.text.y = element_text(size=14))
d
dev.off()

#This isolates the significant SNPs only
CopyBChrom_GWAS_SIG <- subset(CopyBChrom_assoc, neglogp > 7.3)

#I blasted the area surrounding these SNPs and have here a list of all SNPs with greater than 80% identity to the B 
DROP <- c("S2_90555553", "S2_90555560", "S3_5366738", "S6_124809142", "S6_124809191", "S6_124809195", "S8_199282680")

#This drops any SNPs by that name
CopyBChrom_assoc <- CopyBChrom_assoc[-c(which(CopyBChrom_assoc$SNP %in% DROP)),]

pdf("CopyBChrom_GWAS_Filt2.pdf", height=2, width=6)
d <- ggplot(CopyBChrom_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  ggtitle("B Chr. Pseudo Copy Number") +
  geom_hline(yintercept=7.3, color="black", linetype="dotted") +
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  scale_y_continuous(limits=c(0,50), breaks=c(0,10,20,30,40,50)) +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5), legend.position = "none", axis.text.y = element_text(size=14))
d
dev.off()

write.csv(CopyBChrom_assoc, "~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/CopyBChrom_results_PCA.assoc.linear.filt", row.names = FALSE, quote = FALSE)


pdf("All_GWAS.pdf", height=4, width=12)
e <- ggarrange(a,c,b,d)
e
dev.off()

pdf("All_GWAS_Hist.pdf", height=4, width=12)
f <- ggarrange(w,x,y,z)
f
dev.off()


#########################################################################
#This is a control with the Ab10 SNPs 
#########################################################################

Ab10_assoc <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10_AB10_PCA_results.assoc.logistic", sep="")
Ab10_assoc$neglogp <- -log10(Ab10_assoc$P)
Ab10_assoc$Direction <- ifelse(Ab10_assoc$BETA < 0, "Negative", "Positive")
Ab10_assoc <- subset(Ab10_assoc, TEST == "ADD")
#Drop NA
Ab10_assoc <- subset(Ab10_assoc, is.na(BETA) == FALSE)

pdf("Ab10_Ab10_GWAS.pdf", height=2, width=6)
i <- ggplot(Ab10_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  ggtitle(paste("Ab10\n", nrow(Ab10_assoc))) +
  geom_hline(yintercept=7.3, color="black", linetype="dotted") +
  
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5))
i
dev.off()



K10L2_assoc <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/K10L2_K10L2_PCA_results.assoc.logistic", sep="")
K10L2_assoc$neglogp <- -log(K10L2_assoc$P)
K10L2_assoc$Direction <- ifelse(K10L2_assoc$BETA < 0, "Negative", "Positive")
K10L2_assoc <- subset(K10L2_assoc, TEST == "ADD")
#Drop NA
K10L2_assoc <- subset(K10L2_assoc, is.na(BETA) == FALSE)

pdf("K10L2_K10L2_GWAS.pdf", height=2, width=6)
j <- ggplot(K10L2_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  ggtitle(paste("K10L2\n", nrow(K10L2_assoc))) +
  geom_hline(yintercept=7.3, color="black", linetype="dotted") +
  
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5))
j
dev.off()


BChrom_assoc <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/BChrom_B_results.assoc.logistic", sep="")
BChrom_assoc$neglogp <- -log10(BChrom_assoc$P)
BChrom_assoc$Direction <- ifelse(BChrom_assoc$BETA < 0, "Negative", "Positive")
BChrom_assoc <- subset(BChrom_assoc, TEST == "ADD")
#Drop NA
BChrom_assoc <- subset(BChrom_assoc, is.na(BETA) == FALSE)

pdf("Bchr_Bchr_GWAS.pdf", height=2, width=6)
k <- ggplot(BChrom_assoc, aes(x=BP, y=neglogp, color = Direction)) +
  geom_point(alpha=0.5, size=1) +
  facet_wrap(~CHR, ncol=10, nrow=1, scales="free_x") +
  scale_color_manual(values = c("Positive"= "red", "Negative"="blue")) +
  ggtitle(paste("B Chr. \n", nrow(BChrom_assoc))) +
  geom_hline(yintercept=7.3, color="black", linetype="dotted") +
  
  labs(x="Genomic Location", y="-log(p)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), plot.title=element_text(hjust=0.5))
k
dev.off()

pdf("All_GWAS_CDHControl.pdf", height=6, width=6)
l <- ggarrange(i,j,k, ncol=1 )
l
dev.off()
