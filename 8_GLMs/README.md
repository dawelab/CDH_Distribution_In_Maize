
# 8. GLMs
These scripts use generalized linear models (GLM) to explore the relationship between chromosome drive haplotypes and possible causal factors. 

1. This script extracts environmental data
2. This scipt averages the vapor pressue across all 12 months available from worldclim2, an identical script was used to do the same for solar radiation
3. This script checks for correlation between environmental variables
4. This script merges whole genome PCA data from PLINK with the meta data and environmental data. Each CDH has it's own unique whole genome PCA. We maped reads to the Mo17 genome with the relevant CDH appended to reduce mismapping. The csv files are the output. 
5. This script performs the GLMs on Ab10 with and without SNPs and conducts deviance partitionsing. Additional files are the outputs.
6. This script performs the GLMs on B chr with and without SNPs and conducts vdeviance partitionsing. Additional files are the outputs.
7. This script performs the GLMs on K10L2 with and without SNPs and conducts deviance partitionsing. Additional files are the outputs.
8. This script performs the GLMs for B chr pseudo copy number with and without SNPs and conducts deviance partitionsing. Additional files are the outputs.
9. This plots the results of all the GLMs for 5-8.
10. This plots the relationship between all CDHs and B chr pseudo copy number and the elevation at which they were collected.
11. This plots the output of the deviancec partitioning output by 5-8
