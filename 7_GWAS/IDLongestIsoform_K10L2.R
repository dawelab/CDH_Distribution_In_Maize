library(tidyr)
library(dplyr)

#generated .fai in bash using module load SAMtools/1.18-GCC-12.3.0 samtools faidx CI66_K10L2.v2.protein.K10L2Hap.fasta

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder")

K10L2HapProt <- read.delim("CI66_K10L2.v2.protein.K10L2Hap.fasta.fai", header=FALSE)

names(K10L2HapProt) <- c("name", "length", "offset", "linebases", "linewidth")


K10L2HapProt <- separate(K10L2HapProt, col = name, into=c("gene", "isoform"), sep="_")

genes <- unique(K10L2HapProt$gene)

Long <- K10L2HapProt %>% group_by(K10L2HapProt$gene) %>% arrange(-length) %>% slice(1)
length(Long$gene)

Long$ID <- paste(Long$gene, Long$isoform, sep="_")

#Print a list of gene names and numbers of the longest isoforms
write.table(Long$ID, file = "CI66_K10L2.v2.protein.K10L2HapGeneNamesLongest.txt", sep="\n", row.names = FALSE, col.names = FALSE, quote = FALSE)