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
