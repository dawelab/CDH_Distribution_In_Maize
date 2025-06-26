#This file is a filtered version of https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main/Data/B73.all.csv. It contains only genes that are single copy across the pangenome. 
CORE <- read.csv("single_core_B73v5.csv")

write.table(CORE$gene, file="single_core_B73v5_genename.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
