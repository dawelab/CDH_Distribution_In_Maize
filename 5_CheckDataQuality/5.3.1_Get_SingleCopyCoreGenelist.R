CORE <- read.csv("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/single_core_B73v5.csv")

write.table(CORE$gene, file="single_core_B73v5_genename.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
