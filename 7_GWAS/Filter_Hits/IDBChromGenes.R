library(tidyverse)

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder")

#Removed the header to make this file compatible with R 
BChr_GFF <- read.delim("Zm-B73_B_CHROMOSOME-MBSC-1.0_Zm00044a.1.gene.gff3", header = FALSE)
colnames(BChr_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")


#This selects only genes on the BChr haplotype
BChr_GFF <- subset(BChr_GFF, feature == "gene" & seqname == "chrB")

#This separates out the gene ID 
BChr_GFF_GENE_temp1 <- separate(BChr_GFF, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')
BChr_GFF_GENE_temp2 <- separate(BChr_GFF_GENE_temp1, col= ID, into=c("Trash", "ID"), sep="=")
BChr_GFF_GENE_temp3 <- BChr_GFF_GENE_temp2[,c("seqname", "start", "end", "ID")]
colnames(BChr_GFF_GENE_temp3) <- c("BChr_seqname", "BChr_start", "BChr_end", "BChr_ID")
BChr_GFF_GENE_temp3$BChr_ID <-  gsub("gene:", "", BChr_GFF_GENE_temp3$BChr_ID)
BChr_GFF_GENE <- BChr_GFF_GENE_temp3

#This writes the gene names on the BChr haplotype
write.table(BChr_GFF_GENE$BChr_ID, "Zm-B73_B_CHROMOSOME-MBSC-1.0_Zm00044a.1.BChrGeneNames.txt", row.names = FALSE, quote = FALSE)