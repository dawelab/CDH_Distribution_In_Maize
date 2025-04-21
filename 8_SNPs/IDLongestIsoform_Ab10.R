library(tidyr)
library(dplyr)

#generated .fai in bash using module load SAMtools/1.18-GCC-12.3.0 samtools faidx Ab10_HiFi_v2_corrected.gene.v2.protein.Ab10Hap.fa

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder")

Ab10HapProt <- read.delim("Ab10_HiFi_v2_corrected.gene.v2.protein.Ab10Hap.fa.fai", header=FALSE)

names(Ab10HapProt) <- c("name", "length", "offset", "linebases", "linewidth")


Ab10HapProt <- separate(Ab10HapProt, col = name, into=c("gene", "isoform"), sep="_")

genes <- unique(Ab10HapProt$gene)

Long <- Ab10HapProt %>% group_by(Ab10HapProt$gene) %>% arrange(-length) %>% slice(1)
length(Long$gene)

Long$ID <- paste(Long$gene, Long$isoform, sep="_")

#Print a list of gene names and numbers of the longest isoforms
write.table(Long$ID, file = "Ab10_HiFi_v2_corrected.gene.v2.Ab10hapGeneNamesLongest.txt", sep="\n", row.names = FALSE, col.names = FALSE, quote = FALSE)