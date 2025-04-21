library(tidyverse)

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder")


#Removed the header to make this file compatible with R 
B73Ab10_GFF <- read.delim("Ab10_HiFi_v2_corrected.gene.v2.gff3", header = FALSE)
colnames(B73Ab10_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

#141115174 is the location of the R1 gene

#This selects only genes on the Ab10 haplotype
Ab10_GFF <- subset(B73Ab10_GFF, feature == "gene" & seqname == "chr10" & start >= 141115174)

#This separates out the gene ID 
Ab10_GFF_GENE_temp1 <- separate(Ab10_GFF, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')
Ab10_GFF_GENE_temp2 <- separate(Ab10_GFF_GENE_temp1, col= ID, into=c("Trash", "ID"), sep="=")
Ab10_GFF_GENE_temp3 <- Ab10_GFF_GENE_temp2[,c("seqname", "start", "end", "ID")]
colnames(Ab10_GFF_GENE_temp3) <- c("Ab10_seqname", "Ab10_start", "Ab10_end", "Ab10_ID")
Ab10_GFF_GENE_temp3$Ab10_ID <-  gsub("gene:", "", Ab10_GFF_GENE_temp3$Ab10_ID)
Ab10_GFF_GENE <- Ab10_GFF_GENE_temp3

#This writes the gene names on the Ab10 haplotype
write.table(Ab10_GFF_GENE$Ab10_ID, "Ab10_HiFi_v2_corrected.gene.v2.Ab10hapGeneNames.txt", row.names = FALSE, quote = FALSE)