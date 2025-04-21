Y_MAIZE <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Yamakake Table/Yamakake_Appendix_Table_5_R.xlsx")
TEO <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Yamakake Table/Yamakake_Appendix_Table_10_R.xlsx")
M6_MAIZE <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Yamakake Table/McClintock_Appendix_Table6.xlsx")
M9_MAIZE <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Yamakake Table/McClintock_Appendix_Table9.xlsx")
M12_MAIZE <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Yamakake Table/McClintock_Appendix_Table12.xlsx")
M14_MAIZE <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Yamakake Table/McClintock_Appendix_Table14.xlsx")

TeoAb10Freq <- sum(TEO$`10ab-1`+TEO$`10ab-2`)/(sum(TEO$Total_Chromosomes)/2)
TeoK10L2Freq <- sum(TEO$`10L2-m` + TEO$`10L2-s`)/(sum(TEO$Total_Chromosomes)/2)
TeoBFreq <- sum(TEO$`BChromosome-Plants_With`)/(sum(TEO$Total_Chromosomes)/2)


Maize_Ab10 <- (sum(Y_MAIZE$`10ab-1`) + sum(Y_MAIZE$`10ab-2`) + sum(na.omit(M6_MAIZE$Ab10)) + sum(na.omit(M9_MAIZE$Ab10)) + sum(na.omit(M12_MAIZE$Ab10)) + sum(na.omit(M14_MAIZE$Ab10)))

Total_Maize_Ab10 <- (sum(Y_MAIZE$Total_Chromosomes)/2) + (sum(M6_MAIZE$`Total Chr`)/2) + (sum(M9_MAIZE$`Total Chr`)/2) + (sum(M12_MAIZE$`Total Chr`)/2) + (sum(M14_MAIZE$`Number of Plants`))

MaizeAb10Freq <- (Maize_Ab10/Total_Maize_Ab10)*100




Maize_K10L2 <- (sum(Y_MAIZE$`10L2-l`) + sum(Y_MAIZE$`10L2-m`) + sum(Y_MAIZE$`10L2-s`) + sum(na.omit(M6_MAIZE$K10L2 )) + sum(na.omit(M9_MAIZE$`10L2`)) + sum(na.omit(M12_MAIZE$`10L2(m and large)`)))

Total_Maize_K10L2 <- ((sum(Y_MAIZE$Total_Chromosomes)/2) + (sum(M6_MAIZE$`Total Chr`)/2) + (sum(M9_MAIZE$`Total Chr`)/2) + (sum(M12_MAIZE$`Total Chr`)/2))

MaizeK10L2Freq <- (Maize_K10L2/Total_Maize_K10L2)*100





Maize_BChr <- (sum(Y_MAIZE$`BChromosome-Plants_With`) + sum(na.omit(M6_MAIZE$B)) + sum(na.omit(M9_MAIZE$B)) + sum(na.omit(M12_MAIZE$B)) + sum(na.omit(M14_MAIZE$`B Chrom`)))

Total_Maize_BChr <- ((sum(Y_MAIZE$Total_Chromosomes)/2) + (sum(M6_MAIZE$`Total Chr`)/2) + (sum(M9_MAIZE$`Total Chr`)/2) + (sum(M12_MAIZE$`Total Chr`)/2) + (sum(M14_MAIZE$`Number of Plants`)))

MaizeBChrFreq <- (Maize_BChr/Total_Maize_BChr)*100

