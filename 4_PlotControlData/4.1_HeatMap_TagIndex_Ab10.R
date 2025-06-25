library(ggplot2)
library(vroom)
library(tidyverse)
library(data.table)
library(pheatmap)
library(circlize)
library(readxl)
library(dichromat)
library(RColorBrewer)
library(ggpubr)
library(stringr)

#This sets the working directory
setwd("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/4_PlotControlData")

########################################This loads and preps the data 
GROUPS <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

MERGE_Ab10Hap_RPM <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.Ab10.RPM.txt")

#7 columns were classed as logicals because they have all NAs, these are the ones that had 0 coverage

BINS <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Bins_NoOverlap.table.xlsx")

#This drops any column that had an average coverage of 0, resulting in all values being NA
DT <- as.data.table(MERGE_Ab10Hap_RPM)
DT <- DT[,which(unlist(lapply(DT, function(x)!all(is.na(x))))),with=F]
MERGE_Ab10Hap_RPM_FILT_1 <- as.data.frame(DT)

#This saves the list of columns that were dropped
DROPPED <- setdiff(colnames(MERGE_Ab10Hap_RPM), colnames(MERGE_Ab10Hap_RPM_FILT_1))

#This removes any tag with a MAPQ less than 20
MERGE_Ab10Hap_RPM_FILT <- subset(MERGE_Ab10Hap_RPM_FILT_1, MAPQ >= 20)

###This selects only the controls
#This identifies all the index numbers for names marked as controls
col.num <- grep(".DC", colnames(MERGE_Ab10Hap_RPM_FILT))
#This is 188 because 4 controls were dropped as they all had NAs. I believe this is a product of them having naming issues that prevented the correct mean from being pulled in the normalization step. 

#This appends the first 6 columns to that list
col.num <- append(col.num, values = c(1:6))
#This selects only the columns that appear in the above list
MERGE_Ab10Hap_RPM_FILT_CONTROLS <- MERGE_Ab10Hap_RPM_FILT[,sort(c(col.num))]

#This drops the Blanks
DROP <- grep("BLANK", colnames(MERGE_Ab10Hap_RPM_FILT_CONTROLS) )
MERGE_Ab10Hap_RPM_FILT_CONTROLS <- MERGE_Ab10Hap_RPM_FILT_CONTROLS[,-DROP]

#This drops misclassified controls, unknown controls, K10L2 which will be handled in a different model, and B73 which have background issues
DROP_1 <- c("W23_AB10-I.11.DC1", "W23_AB10-I.13.DC1", "W23_AB10-II.36.DC1", "W23_N10.14.DC1")
DROP_2 <- c("B73_N10.1.DC2", "B73_N10.2.DC2", "B73_N10.3.DC2")
K10L2 <- subset(GROUPS, Ab10_Status == "K10L2" | Ab10_Status == "Ab10_Unknown")
DROP_3 <- c(DROP_1, DROP_2, K10L2$Name)
DROP_5 <- which(colnames(MERGE_Ab10Hap_RPM_FILT_CONTROLS)  %in% DROP_3)
MERGE_Ab10Hap_RPM_FILT_CONTROLS <- MERGE_Ab10Hap_RPM_FILT_CONTROLS[,-DROP_5]

#This reduces the number of N10 samples included
GROUPS_CONTROLS_Ab10 <- subset(GROUPS, (GROUPS$Data_Source == "Dawe_Lab_1" | GROUPS$Data_Source == "Dawe_Lab_2") & (GROUPS$Ab10_Status == "Ab10-I" | GROUPS$Ab10_Status == "Ab10-II" | GROUPS$Ab10_Status == "Ab10-III"))

#This selects 5 N10 samples
GROUPS_CONTROLS_N10 <- subset(GROUPS, (GROUPS$Data_Source == "Dawe_Lab_1" | GROUPS$Data_Source == "Dawe_Lab_2") & (GROUPS$B_Chrom_Status == "No" & GROUPS$Ab10_Status == "N10"))

GROUPS_CONTROLS_N10 <- GROUPS_CONTROLS_N10[1:5,]

GROUPS_CONTROLS <- rbind(GROUPS_CONTROLS_Ab10, GROUPS_CONTROLS_N10)

#This identifies all the index numbers for names that appear in the GROUPS control only name fields
col.num <- which(colnames(MERGE_Ab10Hap_RPM_FILT_CONTROLS) %in% GROUPS_CONTROLS$Name)
#This appends the first 6 columns to that list
col.num <- append(col.num, values = c(1:6))
#This selects only the columns that appear in the above list
MERGE_Ab10Hap_RPM_FILT_CONTROLS <- MERGE_Ab10Hap_RPM_FILT_CONTROLS[,sort(c(col.num))]

#this adds a bin value to each line in the MERGE_Ab10Hap_RPM file 
MERGE_Ab10Hap_RPM_FILT_CONTROLS$bin <- NA

i=1
x=1
for (i in 1:nrow(MERGE_Ab10Hap_RPM_FILT_CONTROLS)) {
  print(i)
  for (x in 1:nrow(BINS)) {
    BIN_NUM <- BINS[x,4][[1]]
    START <- BINS[x,2][[1]]
    END <- BINS[x,3][[1]]
    ####CHANGE THIS NUMBER TO 3 WHEN YOU FIX THE PROBLEM WITH THE START
    if (MERGE_Ab10Hap_RPM_FILT_CONTROLS[i,3] > START & MERGE_Ab10Hap_RPM_FILT_CONTROLS[i,3] < END) {
      MERGE_Ab10Hap_RPM_FILT_CONTROLS[i,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS)] <- BIN_NUM
    }
  }
}

#This outputs the bins as well as the tag coordinates to verify correct identification
test <- MERGE_Ab10Hap_RPM_FILT_CONTROLS[,c(1:4,ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS))]

#This function goes over each row and divides each value by the max in that row 
MinMax = function(xx) { sweep(xx, 1, apply(xx, 1, max), '/') }

#This is being a problem and I'm not sure why 
SORT_temp1 <- subset(GROUPS, Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2")
SORT_temp2 <- SORT_temp1[order(SORT_temp1$Ab10_Status),]
SORT_temp3 <- SORT_temp2$Name
ADD <- c("Tag", "Chr", "Start", "End", "MAPQ", "Strand", "bin")
SORT_temp4 <- c(ADD, SORT_temp3)

#This removes the missing samples 
SORT <- SORT_temp4[! SORT_temp4 %in% setdiff(SORT_temp4, colnames(MERGE_Ab10Hap_RPM_FILT_CONTROLS))]

#This manually alters SORT
SORT <- c(SORT[1:7], SORT[c(47,48,49, 8,14,28,35,37,39,42,45)], SORT[c(9:13, 15:27, 29:34, 36,38,40:41,43:44,46,55:59,50:54)], SORT[c(74,83, 60:62,65,68,70,73,76:77,80:81,86:87,91:92)], SORT[c(63:64, 66:67, 69,71:72, 75,78:79, 82,84:85, 88:90, 93:135)])

#This orders the data frame
MERGE_Ab10Hap_RPM_FILT_CONTROLS <- as.data.frame(MERGE_Ab10Hap_RPM_FILT_CONTROLS[,SORT])

MERGE_Ab10Hap_RPM_FILT_CONTROLS_TAGNUM <- MERGE_Ab10Hap_RPM_FILT_CONTROLS

######## Tag Sums
#This sums the number of tag numbers across bins, but excludes the bin column from the summing
A_SUM = as.data.frame(apply(MERGE_Ab10Hap_RPM_FILT_CONTROLS_TAGNUM[,c(8:ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_TAGNUM))], 2, function(xx) { by(xx, MERGE_Ab10Hap_RPM_FILT_CONTROLS_TAGNUM$bin, sum)  }))

######## Tag Density
#This sums the number of tag numbers across bins, but excludes the bin column from the summing
A_DENS = as.data.frame(apply(MERGE_Ab10Hap_RPM_FILT_CONTROLS_TAGNUM[,c(8:ncol(MERGE_Ab10Hap_RPM_FILT_CONTROLS_TAGNUM))], 2, function(xx) { by(xx, MERGE_Ab10Hap_RPM_FILT_CONTROLS_TAGNUM$bin, function(yy) { sum(yy > 0) }) }))

#This takes the square root of each value

B_DENS <- sqrt(A_DENS)

#This combines the density and sum data
B_COMB <- A_SUM+B_DENS

#This scales the data so that the min is always 0 and the max is always 1. 
C <- MinMax(B_COMB)

#This adds bins back
C$bin <- rownames(C)

C_bed <- merge(BINS, C, by="bin")

#This drops the identifier columns to create a numeric matrix
rownames(C_bed) <- C_bed$bin
heat_data2 <- C_bed[,-c(1:4)]

#This generates a data frame with the Ab10 Status of all the lines after sorting them 
Name <- as.data.frame(colnames(heat_data2))
colnames(Name) <- "Name"
df_extra <- merge(Name, GROUPS)
df_col <- df_extra[,c("Name", "Ab10_Status")]
rownames(df_col) = df_col$Name

#This orders the df_col rows to match the data columns
df_col <- as.data.frame(df_col[SORT[-c(1:7)],])
df_col <- as.data.frame(df_col$Ab10_Status)
colnames(df_col) <- c("Chr10 Type")
df_col$`Chr10 Type` <- as.factor(df_col$`Chr10 Type`)


#This reorders the actual heat data to match above 
heat_data2 <- heat_data2[,SORT[-c(1:7)]]

#This line sets the rownames in the df as the column names in the heat map data. IT is very important
rownames(df_col) = colnames(heat_data2)

#This creates a dataframe for the position data 
pos <- C_bed$start
df_row <- as.data.frame(pos)
df_row$feature <- " "
i=1
#This assigns Ab10 regions
for(i in 1:nrow(df_row)) {
  VALUE <- df_row[i,1]
  df_row[i,2] <- ifelse(VALUE >= 142472000 & VALUE <= 146699300, 'TR-1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 148964528 & VALUE <= 149082763, 'Trkin', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 150656000 & VALUE <= 153145000, 'TR-1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 152050000 & VALUE <= 156350000, 'Shared Region', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 158250000 & VALUE <= 166820000, 'Shared Region', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 157485200 & VALUE <= 159356550, 'TR-1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 174433450 & VALUE <= 182846100, 'knob180', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 189326066 & VALUE <= 190330226, 'Kindr', df_row[i,2])
}

df_row <- as.data.frame(df_row$feature)
rownames(df_row) = rownames(C_bed)
colnames(df_row) <- "Feature"
df_row$Feature <- as.factor(df_row$Feature)


colors <- list(`Chr10 Type` = c("Ab10-I" = "#440154FF","Ab10-II" = "#2A788EFF","Ab10-III" =  "#7AD151FF", "N10" ="grey70"), Feature = c('TR-1' = "#56B4E9", 'Trkin' = "#0072B2", 'Shared Region' = "#E69F00", 'knob180' = "#D55E00", 'Kindr' = "#CC79A7", " " = "white" ))

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_Ab10_recolor_publish.pdf", height = 5, width = 8)
pheatmap(heat_data2, scale = "none", show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE, treeheight_row = 0, treeheight_col = 0, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), annotation_col = df_col, main="Ab10 Haplotypes\nTag Index in 1MB bins", fontsize = 10)
dev.off()

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_Ab10_recolor.pdf", height = 5, width = 5)
pheatmap(heat_data2, scale = "none", show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE, treeheight_row = 0, treeheight_col = 0, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), annotation_col = df_col, main="Ab10 Haplotypes\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_Ab10_ColNames_recolor.pdf", height = 5, width = 5)
pheatmap(heat_data2, scale = "none", show_rownames = FALSE, show_colnames = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, treeheight_row = 0, treeheight_col = 0, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), annotation_col = df_col, fontsize_col = 2, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_Ab10_ColClust_recolor.pdf", height = 5, width = 5)
pheatmap(heat_data2, scale = "none", show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = TRUE, clustering_method = "ward.D", treeheight_row = 0, treeheight_col = 50, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), annotation_col = df_col, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_Ab10_ColClust_ColNames_recolor.pdf", height = 5, width = 10)
pheatmap(heat_data2, scale = "none", show_rownames = FALSE, show_colnames = TRUE, cluster_rows = FALSE, cluster_cols = TRUE, clustering_method = "ward.D", treeheight_row = 0, treeheight_col = 50, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), annotation_col = df_col, fontsize_col = 3, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_Ab10_ColClust_RowClust_recolor.pdf", height = 5, width = 10)
pheatmap(heat_data2, scale = "none", show_rownames = FALSE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D", treeheight_row = 50, treeheight_col = 50, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), annotation_col = df_col, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_Ab10_ColClust_RowClust_RowNames_recolor.pdf", height = 5, width = 10)
pheatmap(heat_data2, scale = "none", show_rownames = TRUE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D", treeheight_row = 50, treeheight_col = 50, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), annotation_col = df_col, fontsize_row = 3, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()
