library(vcfR)
library(SNPfiltR)

vcfR <- read.vcfR("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/GWAS_Mo17_K10L2.chr1.vcf.gz")

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/K10L2_VCF_Plots")

popmap<-data.frame(id=colnames(vcfR@gt)[2:length(colnames(vcfR@gt))],pop=substr(colnames(vcfR@gt)[2:length(colnames(vcfR@gt))], 1,15))

#This uses the SNPfiltR package to generate summary plots about the SNPs 
pdf("Depth_Quality.pdf")
hard_filter(vcfR=vcfR)
dev.off()

#This explores missing data
pdf("Missing_By_Sample.pdf")
missing_by_sample(vcfR=vcfR, popmap = popmap)
dev.off()

#This plots missing data by SNP now
pdf("Missing_By_SNP.pdf")
missing_by_snp(vcfR)
dev.off()

#This plots allele balance
pdf("Allele_Balance.pdf")
vcfR<-filter_allele_balance(vcfR)
dev.off()

#This plots the max depth
pdf("MaxDepth.pdf")
max_depth(vcfR)
dev.off()
