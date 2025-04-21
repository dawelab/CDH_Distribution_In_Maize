library(tidyverse)

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder")

#Removed the header to make this file compatible with R 
K10L2_GFF <- read.delim("CI66_K10L2_v1.gene.v2.gff3", header = FALSE)
colnames(K10L2_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

#2736860 is the start of the K10L2 haplotype

#This selects only genes on the K10L2 haplotype
K10L2_GFF <- subset(K10L2_GFF, feature == "gene" & seqname == "K10L2" & start >= 2736860)

#This separates out the gene ID 
K10L2_GFF_GENE_temp1 <- separate(K10L2_GFF, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')
K10L2_GFF_GENE_temp2 <- separate(K10L2_GFF_GENE_temp1, col= ID, into=c("Trash", "ID"), sep="=")
K10L2_GFF_GENE_temp3 <- K10L2_GFF_GENE_temp2[,c("seqname", "start", "end", "ID")]
colnames(K10L2_GFF_GENE_temp3) <- c("K10L2_seqname", "K10L2_start", "K10L2_end", "K10L2_ID")
K10L2_GFF_GENE_temp3$K10L2_ID <-  gsub("gene:", "", K10L2_GFF_GENE_temp3$K10L2_ID)
K10L2_GFF_GENE <- K10L2_GFF_GENE_temp3

#This writes the gene names on the K10L2 haplotype
write.table(K10L2_GFF_GENE$K10L2_ID, "CI66_K10L2_v1.gene.v2.K10L2hapGeneNames.txt", row.names = FALSE, quote = FALSE)