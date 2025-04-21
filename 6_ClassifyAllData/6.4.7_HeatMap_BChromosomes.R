library(ggplot2)
library(vroom)
library(tidyverse)
library(dplyr)
library(data.table)
library(circlize)
library(readxl)
library(dichromat)
library(RColorBrewer)
library(ggpubr)
library(stringr)
library(ComplexHeatmap)

setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/6_ClassifyAllData")

BHigh <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Kmeans_Exp_Groups_High_All.csv")

BLow <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Kmeans_Exp_Groups_Low_All.csv")

#This loads in the filtered, Romero Navarro merged data set. 
B_COMB <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/BChrom/BChrom_TagSumPlusTagDensAllSamples.csv")

#This loads in the groups file with modified names
DF <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/Ab10Model/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table")

#This loads in the BINS file
BINS <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Bins_NoOverlap_BChrom.table.xlsx")

#This loads the original groups file
GROUPS <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")


BHigh_edit <- as.data.frame(apply(BHigh, 2, function(x) {gsub("Yes", 1, x)}))
BHigh_edit <- as.data.frame(apply(BHigh_edit, 2, function(x) {gsub("No", 0, x)}))
BHigh_edit2 <- as.data.frame(apply(BHigh_edit[,-c(1)], 2, function(x) {as.numeric(x)}))
BHigh_edit2$Name <- BHigh_edit$Name
BHigh_edit2$All <- rowSums(BHigh_edit2[,-c(ncol(BHigh_edit2))])
BHigh_edit2$Prop <- BHigh_edit2$All/125
BHigh_edit2$Model <- "High"

ggplot(data=BHigh_edit2, aes(x=Prop)) +
  geom_histogram()
ggsave("B_Chromosome_High_Prop_Calls.png")

ggplot(data=BHigh_edit2, aes(x=Prop)) +
  geom_histogram() +
  ylim(0,400)
ggsave("B_Chromosome_High_Prop_Calls_Zoom.png")

BHigh_final <- subset(BHigh_edit2, Prop >= 0.80)
BHigh_final <- BHigh_final[,c("Name", "Prop", "Model")]




BLow_edit <- as.data.frame(apply(BLow, 2, function(x) {gsub("Yes", 1, x)}))
BLow_edit <- as.data.frame(apply(BLow_edit, 2, function(x) {gsub("No", 0, x)}))
BLow_edit2 <- as.data.frame(apply(BLow_edit[,-c(1)], 2, function(x) {as.numeric(x)}))
BLow_edit2$Name <- BLow_edit$Name
BLow_edit2$All <- rowSums(BLow_edit2[,-c(ncol(BLow_edit2))])
BLow_edit2$Prop <- BLow_edit2$All/125
BLow_edit2$Model <- "Low"

ggplot(data=BLow_edit2, aes(x=Prop)) +
  geom_histogram()
ggsave("B_Chromosome_Low_Prop_Calls.png")

ggplot(data=BLow_edit2, aes(x=Prop)) +
  geom_histogram() +
  ylim(0,400)
ggsave("B_Chromosome_Low_Prop_Calls_Zoom.png")

BLow_final <- subset(BLow_edit2, Prop >= 0.95)
BLow_final <- BLow_final[,c("Name", "Prop", "Model")]

#This brings together all B chromosome positive samples
BChrom_All <- rbind(BHigh_final, BLow_final)

#This writes out the B Chromosome Positive Samples
NoB_final <- subset(BLow_edit2, Prop < 0.95)
NoB_final <- NoB_final[,c("Name", "Prop", "Model")]
NoB_final$Model <- "Neither"

BChrom_All_v2 <- BChrom_All
BChrom_All_v2$KMeans_BChrom <- "Yes"
NoB_final$KMeans_BChrom <- "No"
BChrom_Final <- rbind(BChrom_All_v2, NoB_final)

write.csv(BChrom_Final, file = "BChromosome_Calls_final.csv", row.names = FALSE, quote = FALSE)

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

#This selects only control data, excluding all K10L2, B73, and mis-classified controls. 
DF_EXTRA <- unique(DF[,-c(2:3)])
DF_EXTRA_2 <- subset(DF_EXTRA, (Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2") & Ab10_Status != "K10L2" & Name != "B73_N10.1.DC2" & Name != "B73_N10.2.DC2" & Name != "B73_N10.3.DC2" & Name != "NSL-2833_B-Chrom.2.DC2" & Name != "B542C_L289_B-Chrom.1.DC2")

#This subsets the heat_data file to only the lines selected above
heat_data_control <- heat_data[,DF_EXTRA_2$Name]

#This isolates only the experimental samples
DF_EXTRA_3 <- subset(DF_EXTRA, Data_Source != "Dawe_Lab_1" & Data_Source != "Dawe_Lab_2")

#This selects only B chromosome positive experimental samples
DF_EXTRA_3 <- DF_EXTRA_3[c(which(DF_EXTRA_3$Name %in% BChrom_All$Name)),]

heat_data_B <- heat_data[,DF_EXTRA_3$Name]

#This generates the dataframe to annotate the columns
df_col <- DF_EXTRA[,c("Name", "B_Chrom_Status", "Data_Source")]
colnames(df_col) <- c("Name", "B Chrom", "Data Source")
df_col$`B Chrom` <- as.factor(df_col$`B Chrom`)
df_col$`Data Source` <- as.factor(df_col$`Data Source`)

#This selects df_col for the controls
df_col_control <- df_col[which(df_col$Name %in% DF_EXTRA_2$Name),]
df_col_control$KMeans_B <- df_col_control$`B Chrom`

#This selects df_col for the Bs
df_col_B <- df_col[which(df_col$Name %in% DF_EXTRA_3$Name),]
df_col_B$KMeans_B <- "Yes"

#This brings together both df_col
df_col_all <- rbind(df_col_control, df_col_B)
df_col_all$`B Chrom` <- factor(df_col_all$`B Chrom`, levels=c("Yes", "No", "Unknown"))

#This rbrings together the heat data
heat_data_all <- cbind(heat_data_control, heat_data_B)
heat_data_all <- heat_data_all[,c(df_col_all$Name)]

#This generates the dataframe to annotate the rows
pos <- C_bed$start
df_row <- as.data.frame(pos)
df_row$feature <- " "

#This assigns Ab10 regions
for(i in 1:nrow(df_row)) {
  VALUE <- df_row[i,1]
  #df_row[i,2] <- ifelse(VALUE >= 1 & VALUE <= 237770, 'B Short Arm', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 237771 & VALUE <= 1113960, 'B Repeat', df_row[i,2])
  #df_row[i,2] <- ifelse(VALUE >= 1142820 & VALUE <= 8162450, 'Proximal Heterochromatin', df_row[i,2])
  #df_row[i,2] <- ifelse(VALUE >= 8162451 & VALUE <= 20739110, 'Proximal euchromatin 1', df_row[i,2])
  #df_row[i,2] <- ifelse(VALUE >= 20739111 & VALUE <= 43049286, 'Proximal euchromatin 2', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 1113960 & VALUE <= 2866338, 'knob180', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 43049287 & VALUE <= 44287294, 'Distal Heterochromatin', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 44287295 & VALUE <= 63186011, 'Distal Heterochromatin', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 63186012 & VALUE <= 71189000, 'Distal Heterochromatin', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 71189000 & VALUE <= 101843150, 'Distal Heterochromatin', df_row[i,2])
  #df_row[i,2] <- ifelse(VALUE >= 101843151 & VALUE <= 106643106, 'distal euchromatin', df_row[i,2])
}

df_row <- as.data.frame(df_row$feature)
rownames(df_row) = rownames(C_bed)
colnames(df_row) <- "Feature"
df_row$Feature <- factor(df_row$Feature, levels= c('B Repeat','knob180', 'Distal Heterochromatin'))


#This determines the column annotation colors 
columnBStat_colors <- c("Yes" = "black", "No" ="grey70", "Unknown" = "white")

#columndatasource_colors <- c("Dawe_Lab_1" = "gold", "Dawe_Lab_2" = "limegreen", "Romero-Navarro_etal_2017" = "#008080", "Swarts_etal_2017" = "#310062", "Romay_etal_2013" = "mediumorchid", "NA" = "white")

row_colors <- c('B Repeat' = "lightpink2", 'knob180' = "#D55E00",'Distal Heterochromatin' =  "#009e73")

#This determines the grid colors
grid_color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500)

#This makes the heatmap annotations
ha_col = HeatmapAnnotation(`B Chrom` = df_col_all$`B Chrom`,
                           col = list(`B Chrom` = columnBStat_colors))

ha_row = rowAnnotation(Feature = df_row$Feature,
                       col = list(Feature = row_colors),
                       na_col = "white")

#Thsi reorders the samples
df_col_all$KMeans_B <- factor(df_col_all$KMeans_B, levels=c("Yes", "No"))


#This rescales the data only for plotting, it does not affect grouping. Without this the differences are hard to see 

pdf(file = "HeatMap_All_BChrom.pdf", height = 8, width = 8)
a <- Heatmap(as.matrix(heat_data_all),
             #plotting
             col = grid_color,
             column_title = "Tag Index Across B Chromosome",
             column_title_gp = gpar(fontsize = 20, fontface = "bold"),
             row_title = "B Chromosome",
             border_gp = gpar(col = "black", lty = 1),
             
             #clustering
             cluster_rows = FALSE,
             cluster_columns = TRUE,
             column_split = df_col_all$KMeans_B,
             clustering_method_columns = "ward.D",
             column_gap = unit(5, "mm"),
             
             #annotations
             top_annotation = ha_col,
             left_annotation = ha_row,
             show_row_names = FALSE,
             show_column_names = FALSE,
             column_names_gp = grid::gpar(fontsize = 1),
)
print(a)
dev.off()


#This rescales the data
heat_data_all_rescale <- MinMax(heat_data_all)

pdf(file = "HeatMap_All_BChrom_rescale.pdf", height = 10, width = 13)
a <- Heatmap(as.matrix(heat_data_all_rescale),
             #plotting
             col = grid_color,
             column_title = "Tag Index Across B Chromosome",
             column_title_gp = gpar(fontsize = 20, fontface = "bold"),
             row_title = "B Chromosome",
             border_gp = gpar(col = "black", lty = 1),
             
             #clustering
             cluster_rows = FALSE,
             cluster_columns = TRUE,
             column_split = df_col_all$KMeans_B,
             clustering_method_columns = "ward.D",
             column_gap = unit(5, "mm"),
             
             #annotations
             top_annotation = ha_col,
             left_annotation = ha_row,
             show_row_names = FALSE,
             show_column_names = TRUE,
             column_names_gp = grid::gpar(fontsize = 1),
)
print(a)
dev.off()

#This rescales again
heat_data_all_rescale2 <- log(heat_data_all+1)

pdf(file = "HeatMap_All_BChrom_rescale2.pdf", height = 8, width = 8)
a <- Heatmap(as.matrix(heat_data_all_rescale2),
             #plotting
             col = grid_color,
             column_title = "Tag Index Across B Chromosome",
             column_title_gp = gpar(fontsize = 20, fontface = "bold"),
             row_title = "B Chromosome",
             border_gp = gpar(col = "black", lty = 1),
             
             #clustering
             cluster_rows = FALSE,
             cluster_columns = TRUE,
             column_split = df_col_all$KMeans_B,
             clustering_method_columns = "ward.D",
             column_gap = unit(5, "mm"),
             
             #annotations
             top_annotation = ha_col,
             left_annotation = ha_row,
             show_row_names = FALSE,
             show_column_names = FALSE,
             column_names_gp = grid::gpar(fontsize = 1),
)
print(a)
dev.off()


#This selects only the B chromosomes
df_col_Bonly <- subset(df_col_all, KMeans_B == "Yes")
heat_data_Bonly <- heat_data_all[,which(colnames(heat_data_all) %in% df_col_Bonly$Name)]

#This makes the heatmap annotations
ha_col = HeatmapAnnotation(`B Chrom` = df_col_Bonly$`B Chrom`,
                           `Data Source` = df_col_Bonly$`Data Source`,
                           col = list(`B Chrom` = columnBStat_colors, `Data Source` = columndatasource_colors))

ha_row = rowAnnotation(Feature = df_row$Feature,
                       col = list(Feature = row_colors),
                       na_col = "white")

#This rescales the data
heat_data_Bonly_rescale <- MinMax(heat_data_Bonly)
#This determines the grid colors
grid_color=rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrBr")))(500))

pdf(file = "HeatMap_All_BChrom_BOnly_rescale.pdf", height = 10, width = 13)
a <- Heatmap(as.matrix(heat_data_Bonly_rescale),
             #plotting
             col = grid_color,
             column_title = "Tag Index Across B Chromosome",
             column_title_gp = gpar(fontsize = 20, fontface = "bold"),
             row_title = "B Chromosome",
             border_gp = gpar(col = "black", lty = 1),
             
             #clustering
             cluster_rows = FALSE,
             cluster_columns = TRUE,
             column_split = df_col_Bonly$KMeans_B,
             clustering_method_columns = "ward.D",
             column_gap = unit(5, "mm"),
             
             #annotations
             top_annotation = ha_col,
             left_annotation = ha_row,
             show_row_names = FALSE,
             show_column_names = TRUE,
             column_names_gp = grid::gpar(fontsize = 1),
)
print(a)
dev.off()

pdf(file = "HeatMap_All_BChrom_BOnly_rescale_rowclust.pdf", height = 10, width = 13)
a <- Heatmap(as.matrix(heat_data_Bonly_rescale),
             #plotting
             col = grid_color,
             column_title = "Tag Index Across B Chromosome",
             column_title_gp = gpar(fontsize = 20, fontface = "bold"),
             row_title = "B Chromosome",
             border_gp = gpar(col = "black", lty = 1),
             
             #clustering
             cluster_rows = TRUE,
             cluster_columns = TRUE,
             column_split = df_col_Bonly$KMeans_B,
             clustering_method_columns = "ward.D",
             column_gap = unit(5, "mm"),
             
             #annotations
             top_annotation = ha_col,
             left_annotation = ha_row,
             show_row_names = TRUE,
             show_column_names = TRUE,
             column_names_gp = grid::gpar(fontsize = 1),
             row_names_gp = grid::gpar(fontsize = 5),
)
print(a)
dev.off()





forPCA <- t(heat_data_Bonly_rescale)

#This performs the pca
dist <- dist(forPCA)
pca <- cmdscale(dist, k=10, eig = TRUE)

#This calculates the percent of the variation they account for 
TOT <- sum(pca$eig)
pca$perc_var <- (pca$eig/TOT)*100 

#This plots the eigen vector values. This tells you how much variability is accounted for by each PCA
pdf("PCA_ScreePlot_BChromOnly.pdf", height = 5, width = 5)
plot(pca$perc_var)
dev.off()

df_out <- as.data.frame(matrix(nrow = length(rownames(forPCA)), ncol = 8))
names(df_out) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "Name")
df_out$PC1 <- pca$points[,1]
df_out$PC2 <- pca$points[,2]
df_out$PC3 <- pca$points[,3]
df_out$PC4 <- pca$points[,4]
df_out$PC5 <- pca$points[,5]
df_out$PC6 <- pca$points[,6]
df_out$PC7 <- pca$points[,7]
df_out$Name <-  c(as.character(row.names(forPCA)))

df_out <- merge(df_out, BChrom_All, by = "Name", all.x = TRUE)

#This merges the file with the group information and the PC information
df_out <- merge(df_out, GROUPS, by = "Name")

df_out_controls <- subset(df_out, Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2")
df_out_controls$Control <- "Yes"
df_out_controls <- df_out_controls[,c("Name", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "Data_Source", "Accession", "Group", "Maize_Type", "Ab10_Status", "B_Chrom_Status", "Latitude", "Longitude", "Altitude", "Control", "Model")]

df_out_exp <- subset(df_out, Data_Source != "Dawe_Lab_1" | Data_Source != "Dawe_Lab_2")
df_out_exp$Control <- "No"
df_out_exp <- df_out_exp[,c("Name", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "Data_Source", "Accession", "Group", "Maize_Type", "Ab10_Status", "B_Chrom_Status", "Latitude", "Longitude", "Altitude", "Control", "Model")]

df_out_edit <- rbind(df_out_controls, df_out_exp)

BChrom_colors <- c("blue", "red", "grey")
data_colors <- c("gold", "limegreen", "mediumorchid", "#008080", "#310062")
shapes <- c(8, 15, 16, 3, 4, 7, 17, 18)

#This loop plots the first two PCAs against each other
pdf("PCA1_v_PCA2_BChrom.pdf", height = 8, width = 15)
ggplot() +
  geom_point(data=df_out_edit, alpha = 0.5, aes(x=PC1, y=PC2, size= Control, color = Model, shape = Maize_Type)) +
  ggtitle("PCA of GBS Tags\n on B Chromosome\n B Chromosome Positive Only") +
  scale_shape_manual(values= shapes) +
  scale_color_manual(values = BChrom_colors) +
  labs(color="Model", shape = "Maize Type") +
  #labs(x="PC1 32%", y="PC2 18%") +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()

#This loop plots the first two PCAs against each other
pdf("PCA1_v_PCA2_BChrom_data.pdf", height = 8, width = 15)
ggplot() +
  geom_point(data=df_out_edit, alpha = 0.5, aes(x=PC1, y=PC2, color=Data_Source , shape = Maize_Type)) +
  ggtitle("PCA of GBS Tags\n on B Chromosome\n B Chromosome Positive Only") +
  scale_shape_manual(values= shapes) +
  scale_color_manual(values = data_colors) +
  labs(color="Data Source", shape = "Maize Type") +
  #labs(x="PC1 32%", y="PC2 18%") +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()

pdf("PCA1_v_PCA3_BChrom.pdf", height = 8, width = 15)
ggplot() +
  geom_point(data=df_out_edit, alpha = 0.5, aes(x=PC1, y=PC3, size= Control, color = Model, shape = Maize_Type)) +
  ggtitle("PCA of GBS Tags\n on B Chromosome\n B Chromosome Positive Only") +
  scale_shape_manual(values= shapes) +
  scale_color_manual(values = BChrom_colors) +
  labs(color="Model", shape = "Maize Type") +
  #labs(x="PC1 32%", y="PC2 18%") +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()

pdf("PCA1_v_PCA3_BChrom_data.pdf", height = 8, width = 15)
ggplot() +
  geom_point(data=df_out_edit, alpha = 0.5, aes(x=PC1, y=PC3, color=Data_Source , shape = Maize_Type)) +
  ggtitle("PCA of GBS Tags\n on B Chromosome\n B Chromosome Positive Only") +
  scale_shape_manual(values= shapes) +
  scale_color_manual(values = data_colors) +
  labs(color="Control", shape = "Maize Type") +
  #labs(x="PC1 32%", y="PC2 18%") +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()

pdf("PCA1_v_PCA4_BChrom.pdf", height = 8, width = 15)
ggplot() +
  geom_point(data=df_out_edit, alpha = 0.5, aes(x=PC1, y=PC4, size= Control, color = Model, shape = Maize_Type)) +
  ggtitle("PCA of GBS Tags\n on B Chromosome\n B Chromosome Positive Only") +
  scale_shape_manual(values= shapes) +
  scale_color_manual(values = BChrom_colors) +
  labs(color="Model", shape = "Maize Type") +
  #labs(x="PC1 32%", y="PC2 18%") +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()

pdf("PCA1_v_PCA4_BChrom_data.pdf", height = 8, width = 15)
ggplot() +
  geom_point(data=df_out_edit, alpha = 0.5, aes(x=PC1, y=PC4, color=Data_Source , shape = Maize_Type)) +
  ggtitle("PCA of GBS Tags\n on B Chromosome\n B Chromosome Positive Only") +
  scale_shape_manual(values= shapes) +
  scale_color_manual(values = data_colors) +
  labs(color="Control", shape = "Maize Type") +
  #labs(x="PC1 32%", y="PC2 18%") +
  theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
dev.off()


