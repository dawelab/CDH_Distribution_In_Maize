# 7 GWAS
These scripts call filter, and impute SNPs before running a GWAS for all CDHs and B Chr. pseudo copy number, filtering the resulting hits and investigating the remaining significantly assocaited loci . 

## Call SNPs
1. This script edits the key file for SNP calling. For the SNP calling we merged technical replicates (from Romero-Navarro) immediately as the data is being read in. The final key used is included here.
2. This script builds the genomes used here. We chose to append each CDH one at a time to the Mo17 genome https://doi.org/10.1038/s41588-023-01419-6 to avoid mismapping. We truncated Mo17 chr10 at the colored 1 gene Zm00001eb429330. The Ab10 and K10L2 genomes refered to are from  https://doi.org/10.1093/genetics/iyaf091. The B chrom assembly is from https://doi.org/10.1073/pnas.2104254118.
3. These scripts run TASSEL and call SNPs for each CDH using the edited key and their respective genomes.
4. These scripts isolate landrace samples only and filter the SNPs to remove low quality sites.
5. These scripts split the filtered VCF file indo chunks for faster processing.
6. These scripts impute SNPs on Mo17 chr1-10 based on the population using BEAGLE.
7. These scripts impute SNPs on the CDH itself based on the population using BEAGLE. This was done as a control


## Run GWAS
1. Run GWAS using Plink using SNPs on chr 1-10
2. Run GWAS using Plink using SNPs on the CDH itself as a control


## Filter Hits
Remove_CDH_Orthologs. These scripts are explained in doi: 10.1002/pld3.567 and serve to identify genes in chromosomes 1-10 of the genome that are orthologs to genes on each CDH and may cause mismapping and incorrect associations. 
1. This script identifies SNPs that overlap with transposable elements. The output is used to filter them. 
2. This script checks if the significantly associated SNPs are in genes that have orthologs on the CDH. 

## Explore Hits
1. This script plots the output of the GWAS and performs some additional filters using the output from Filter hits.
2. This script wass used to extract the surrounding sequence from a SNP to use web based BLAST to investigate it
3. This script determined which annotated Mo17 genes the final associated SNPs overlapped with. 
