library(tidyr)
library(dplyr)

#generated .fai in bash using module load SAMtools/1.18-GCC-12.3.0 samtools faidx Zm-B73_B_CHROMOSOME-MBSC-1.0_Zm00044a.1.protein.BChr.fa

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder")

BChrProt <- read.delim("Zm-B73_B_CHROMOSOME-MBSC-1.0_Zm00044a.1.protein.BChr.fa.fai", header=FALSE)

names(BChrProt) <- c("name", "length", "offset", "linebases", "linewidth")


BChrProt <- separate(BChrProt, col = name, into=c("gene", "isoform"), sep="_")

genes <- unique(BChrProt$gene)

Long <- BChrProt %>% group_by(BChrProt$gene) %>% arrange(-length) %>% slice(1)
length(Long$gene)

Long$ID <- paste(Long$gene, Long$isoform, sep="_")

#Print a list of gene names and numbers of the longest isoforms
write.table(Long$ID, file = "Zm-B73_B_CHROMOSOME-MBSC-1.0_Zm00044a.1.protein.BChrGeneNamesLongest.txt", sep="\n", row.names = FALSE, col.names = FALSE, quote = FALSE)