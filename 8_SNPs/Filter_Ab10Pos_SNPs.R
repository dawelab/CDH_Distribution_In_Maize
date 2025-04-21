library(vcfR)
library(SNPfiltR)

vcfR <- read.vcfR("/Volumes/Transcend/AllData_v5.Ab10Hap.AllAb10Pos.vcf.gz")

popmap<-data.frame(id=colnames(vcfR@gt)[2:length(colnames(vcfR@gt))],pop=substr(colnames(vcfR@gt)[2:length(colnames(vcfR@gt))], 1,15))

#This uses the SNPfiltR package to generate summary plots about the SNPs 
hard_filter(vcfR=vcfR)

#This applys a filer 
vcfR<-hard_filter(vcfR=vcfR, depth = 3, gq = 60)

#This filters heterozygous SNPs if they are not roughly equal in frequency
vcfR<-filter_allele_balance(vcfR)

#This plots the max depth
max_depth(vcfR)

#This filters anything with a max depth more than 40 reads
vcfR<-max_depth(vcfR, maxdepth = 20)

vcfR
#After these hard filters I have 7514 SNPs remaining

#This explores missing data
missing_by_sample(vcfR=vcfR, popmap = popmap)

#This sets a filter to remove any sample with more than 90 % missing data
vcfR<-missing_by_sample(vcfR=vcfR, cutoff = .9)

#This drops any SNP that is no longer polymorphic without dropped individuals
popmap<-popmap[popmap$id %in% colnames(vcfR@gt),]
#In order to do this I need to extract only biallelic sites
vcfR <- filter_biallelic(vcfR)

vcfR<-min_mac(vcfR, min.mac = 1)

#This checks to ensure that missing data is not driving clustering patterns 
miss<-assess_missing_data_pca(vcfR=vcfR, popmap = popmap, thresholds = .8, clustering = FALSE)

#This plots missing data by SNP now
missing_by_snp(vcfR)

miss<-assess_missing_data_pca(vcfR=vcfR, popmap = popmap, thresholds = c(0.6 ,.7,.75,.8), clustering = FALSE)

#This filters the data to remove SNPs with more than 50% missing data
vcfR<-missing_by_snp(vcfR, cutoff = .6)

#This filters to remove any SNP with more than two minor alleles
vcfR.mac<-min_mac(vcfR = vcfR, min.mac = 2)

#Looks at clustering with and without minor allele filtering
miss<-assess_missing_data_tsne(vcfR, popmap, clustering = FALSE)
miss<-assess_missing_data_tsne(vcfR.mac, popmap, clustering = FALSE)

#Plots the depth 
dp <- extract.gt(vcfR, element = "DP", as.numeric=TRUE)
heatmap.bp(dp, rlabels = FALSE)

#Plots genotype quality
gq <- extract.gt(vcfR, element = "GQ", as.numeric=TRUE)
heatmap.bp(gq, rlabels = FALSE)

#writes out the final fintered VCF file
vcfR::write.vcf(vcfR, "/Volumes/Transcend/AllData_v5.Ab10Hap.AllAb10Pos.filt.vcf.gz")
