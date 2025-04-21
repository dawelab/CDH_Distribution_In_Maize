
VAR <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/PCA_PercVar_longer.csv")
colnames(VAR) <- c("PC", "Percent Var")

VAR_SUB <- VAR[c(1:10),]
sum(VAR_SUB$`Percent Var`)
#11.35258