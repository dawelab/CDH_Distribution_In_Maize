#module load R/4.4.1-foss-2022b

library(tidyverse)
library(reshape2)
library(splitstackshape)
library(readxl)
library(ggplot2)
library(pafr)
library(Rsamtools)

ORTHO <- read.delim("/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/Mo17_Ab10_proteomes/OrthoFinder/Results_Feb20/Orthologues/Orthologues_Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.protein.Longest/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.protein.Longest__v__Ab10_HiFi_v2_corrected.gene.v2.protein.Ab10Hap.longest.tsv")

Mo17_GFF <- read.delim("/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.nohead.gff3", header = FALSE)
colnames(Mo17_GFF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

#This section alters the Mo17 GFF to be compatible with the OrthoFinder Output
Mo17_GFF_GENE <- subset(Mo17_GFF, feature == "gene")

Mo17_GFF_GENE_temp1 <- separate(Mo17_GFF_GENE, col=attribute, into= c("ID", "biotype", "logic_name"), sep= ';')

Mo17_GFF_GENE_temp2 <- separate(Mo17_GFF_GENE_temp1, col= ID, into=c("Trash", "ID"), sep="=")

Mo17_GFF_GENE_temp3 <- Mo17_GFF_GENE_temp2[,c("seqname", "start", "end", "ID")]
colnames(Mo17_GFF_GENE_temp3) <- c("Mo17_seqname", "Mo17_start", "Mo17_end", "Mo17_ID")

Mo17_GFF_GENE <- Mo17_GFF_GENE_temp3

#This section converts the file so that there is only one gene in each file 
ORTHO_MELT_temp1 <- cSplit(ORTHO, "Zm.Mo17.REFERENCE.CAU.2.0_Zm00014ba.protein.Longest", sep = ",", direction = "long")

ORTHO_MELT_temp2 <- cSplit(ORTHO_MELT_temp1, "Ab10_HiFi_v2_corrected.gene.v2.protein.Ab10Hap.longest", sep = ",", direction = "long")

ORTHO_MELT_temp3 <- separate(ORTHO_MELT_temp2, col = Zm.Mo17.REFERENCE.CAU.2.0_Zm00014ba.protein.Longest, into = c("Mo17_ID", "Mo17_Iso"))

ORTHO_MELT_temp4 <- separate(ORTHO_MELT_temp3, col = Ab10_HiFi_v2_corrected.gene.v2.protein.Ab10Hap.longest , into = c("Ab10_ID", "Ab10_Iso"))

ORTHO_MELT <- ORTHO_MELT_temp4


#This section merges the GFF files with the OrthoFinder Files
DATA <- merge(ORTHO_MELT, Mo17_GFF_GENE)

#This writes out the file 
write.table(DATA, "/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/Mo17_Ab10_proteomes/Mo17GenesOrthologousToAb10.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)