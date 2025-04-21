library(ggplot2)
library(vroom)
library(tidyverse)
library(data.table)
library(circlize)
library(readxl)
library(dichromat)
library(RColorBrewer)
library(ggpubr)
library(stringr)
library(ComplexHeatmap)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/6_ClassifyAllData")

GROUPS <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2Complete.csv")

#This loads in the K10L2 tag index
B_COMB <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/Ab10Model/Ab10_TagSumPlusTagDensAllSamples.csv")

#This loads in the groups file with modified names
DF <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/Ab10Model/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table")

#This loads in the BINS file
BINS <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Bins_NoOverlap.table.xlsx")

#This function goes over each row and divides each value by the max in that row 
MinMax = function(xx) { sweep(xx, 1, apply(xx, 1, max), '/') }

#This scales the data so that the min is always 0 and the max is always 1. 
C <- MinMax(B_COMB)

#This adds bins back
C$bin <- rownames(C)

C_bed <- merge(BINS, C, by="bin")

#This drops the identifier columns to create a numeric matrix
rownames(C_bed) <- C_bed$bin
heat_data <- C_bed[,-c(1:4)]

#This selects only unique entries, columns 2 and 3 include the extra Romero Navarro names to differentiate each replicate and need to be removed to find unique lines. 
DF_EXTRA <- unique(DF[,-c(2:3)])
#This selects only controls and drops samples that are known to be misclassified
DF_EXTRA_2 <- subset(DF_EXTRA, (Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2") & Ab10_Status != "Ab10_Unknown" & Ab10_Status != "K10L2" & Name != "W23_AB10-I.11.DC1"& Name != "W23_AB10-I.13.DC1" & Name != "W23_AB10-II.36.DC1" & Name != "W23_AB10-II.37.DC1" & Name != "W23_N10.14.DC1")

#This subsets only the controls 
heat_data_control <- heat_data[,DF_EXTRA_2$Name]

#Select only the relevant columns
df_col_control <- DF_EXTRA_2[,c("Name", "Ab10_Status", "Data_Source")]
colnames(df_col_control) <- c("Name", "Chr10 Type", "Data Source")
df_col_control$`Chr10 Type` <- as.factor(df_col_control$`Chr10 Type`)
df_col_control$`Data Source` <- as.factor(df_col_control$`Data Source`)

#This adds a K10L2 status
df_col_Ab10 <- subset(df_col_control, `Chr10 Type` != "N10" & `Chr10 Type` != "K10L2")
df_col_Ab10$Ab10_Status_G <- "Ab10"

df_col_N10 <- subset(df_col_control, `Chr10 Type` == "N10")
df_col_N10$Ab10_Status_G <- "N10"

#This selects the Experimental lines
GROUPS_Ab10 <- subset(GROUPS, KMeans_Ab10 == "Ab10")
df_col_exp <- GROUPS_Ab10[,c("Name", "Ab10_Status", "Data_Source")]
colnames(df_col_exp) <- c("Name", "Chr10 Type", "Data Source")
df_col_exp$Ab10_Status_G <- "Ab10"

df_col_all <- rbind(df_col_Ab10, df_col_N10, df_col_exp)

#This section reorders the data by their chromosome 10 status
df_col_all <- df_col_all[order(df_col_all$`Chr10 Type`),]

#This line reorders the heat_data to match the df_col
heat_data <- heat_data[,df_col_all$Name]


#This creates a dataframe for the position data 
pos <- C_bed$start
df_row <- as.data.frame(pos)
df_row$feature <- " "

#This assigns Ab10 regions
for(i in 1:nrow(df_row)) {
  VALUE <- df_row[i,1]
  df_row[i,2] <- ifelse(VALUE >= 142472000 & VALUE <= 146699300, 'TR-1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 148964528 & VALUE <= 149082763, 'Trkin', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 150656000 & VALUE <= 153145000, 'TR-1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 152050000 & VALUE <= 156350000, 'Shared Region', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 158250000 & VALUE <= 166820000, 'Shared Region', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 157485200 & VALUE <= 159356550, 'TR-1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 174433450 & VALUE <= 182846100, 'knob 180', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 189326066 & VALUE <= 190330226, 'Kindr', df_row[i,2])
}

df_row <- as.data.frame(df_row$feature)
rownames(df_row) = rownames(C_bed)
colnames(df_row) <- "Feature"
df_row$Feature <- as.factor(df_row$Feature)

#This determines the column annotation colors 
columntype_colors <- c("Ab10-I" = "#440154FF","Ab10-II" = "#2A788EFF","Ab10-III" =  "#7AD151FF", "N10" ="grey70", "Unknown" = "white")

#columndatasource_colors <- c("Dawe_Lab_1" = "gold", "Dawe_Lab_2" = "limegreen", "Romero-Navarro_etal_2017" = "#008080", "Swarts_etal_2017" = "#310062", "Romay_etal_2013" = "mediumorchid", "NA" = "white")

row_colors <- c('TR-1' = "#56B4E9", 'Trkin' = "#0072B2", 'Shared Region' = "#E69F00", 'knob 180' = "#D55E00", 'Kindr' = "#CC79A7")

#This determines the grid colors
grid_color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500)


#This makes the heatmap annotations and plots it
ha_col = HeatmapAnnotation(`Chr10 Type` = df_col_all$`Chr10 Type`,
                           col = list(`Chr10 Type` = columntype_colors))

ha_row = rowAnnotation(Feature = df_row$Feature,
                       col = list(Feature = row_colors),
                       na_col = "white"
)

df_col_all$Ab10_Status_G <- factor(df_col_all$Ab10_Status_G, levels=c("N10", "Ab10"))

pdf(file = paste("Ab10_Controls_AllPositiveExp", "pdf", sep="."), height = 8, width = 8)
a <- Heatmap(as.matrix(heat_data),
             #plotting
             col = grid_color,
             column_title = "Tag Index Across Ab10 Haplotype",
             column_title_gp = gpar(fontsize = 20, fontface = "bold"),
             row_title = "Ab10 Haplotype", 
             border_gp = gpar(col = "black", lty = 1),
             
             #clustering
             cluster_rows = FALSE, 
             cluster_columns = TRUE,
             column_split = df_col_all$Ab10_Status_G,
             clustering_method_columns = "ward.D",
             column_gap = unit(5, "mm"),
             
             #annotations
             top_annotation = ha_col,
             left_annotation = ha_row, 
             show_row_names=FALSE,
             show_column_names = FALSE,
             column_names_gp = grid::gpar(fontsize = 1),
)
print(a)
dev.off()


Inbred <- heat_data[,c(1:168,246:247)]
df_col_all <- df_col_all[c(1:168, 246:247),]

df_col_all[4,246] <- "Ab10"
df_col_all[4,247] <- "Ab10"
#This makes the heatmap annotations and plots it
ha_col = HeatmapAnnotation(`Chr10 Type` = df_col_all$`Chr10 Type`,
                           col = list(`Chr10 Type` = columntype_colors))


pdf(file = paste("Ab10_Controls_InbredPositive", "pdf", sep="."), height = 8, width = 16)
a <- Heatmap(as.matrix(Inbred),
             #plotting
             col = grid_color,
             column_title = "Tag Index Across Ab10 Haplotype",
             column_title_gp = gpar(fontsize = 20, fontface = "bold"),
             row_title = "Ab10 Haplotype", 
             border_gp = gpar(col = "black", lty = 1),
             
             #clustering
             cluster_rows = FALSE, 
             cluster_columns = TRUE,
             column_split = df_col_all$Ab10_Status_G,
             clustering_method_columns = "ward.D",
             column_gap = unit(5, "mm"),
             
             #annotations
             top_annotation = ha_col,
             left_annotation = ha_row, 
             show_row_names=TRUE,
             show_column_names = TRUE,
             column_names_gp = grid::gpar(fontsize = 5),
)
print(a)
dev.off()
