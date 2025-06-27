#!/usr/bin/Rscript
#module load R/4.4.1-foss-2022b
library(vcfR)
library(devtools)
devtools::install_github("DevonDeRaad/SNPfiltR")
library(SNPfiltR)

setwd("/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/")

vcfR <- read.vcfR("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/CoordinateOnly_Ab10BChrom.chr1to9.filt2.vcf.gz")

popmap<-data.frame(id=colnames(vcfR@gt)[2:length(colnames(vcfR@gt))],pop=substr(colnames(vcfR@gt)[2:length(colnames(vcfR@gt))], 1,15))

#This filters the data to remove SNPs with more than 75% missing data
vcfR<-missing_by_snp(vcfR, cutoff = 0.75)

#This sets a filter to remove any sample with more than 90 % missing data
vcfR<-missing_by_sample(vcfR=vcfR, cutoff = 0.9)

#This filters heterozygous SNPs if they are not roughly equal in frequency
vcfR<-filter_allele_balance(vcfR)

#writes out the final fintered VCF file
vcfR::write.vcf(vcfR, "/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/CoordinateOnly_Ab10BChrom.chr1to9.filt4.vcf.gz")
