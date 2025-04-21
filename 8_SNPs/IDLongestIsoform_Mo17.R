module load R/4.4.1-foss-2022b
R 

library(tidyr)
library(dplyr)

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder")

#This file is from https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-2.0/
Mo17Prot <- read.delim("Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.protein.fa.gz.fai", header=FALSE)

names(Mo17Prot) <- c("name", "length", "offset", "linebases", "linewidth")

Mo17Prot <- separate(Mo17Prot, col = name, into=c("gene", "isoform"), sep="_")

genes <- unique(Mo17Prot$gene)

Long <- Mo17Prot %>% group_by(Mo17Prot$gene) %>% arrange(-length) %>% slice(1)
length(Long$gene)

Long$ID <- paste(Long$gene, Long$isoform, sep="_")

#Print a list of gene names and numbers of the longest isoforms
write.table(Long$ID, file = "Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.protein.GeneNamesLongest.txt", sep="\n", row.names = FALSE, col.names = FALSE, quote = FALSE)