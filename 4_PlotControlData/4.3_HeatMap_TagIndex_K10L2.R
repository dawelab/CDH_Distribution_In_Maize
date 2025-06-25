library(ggplot2)
library(vroom)
library(tidyverse)
library(data.table)
library(pheatmap)
library(readxl)
library(dichromat)
library(RColorBrewer)
library(ggpubr)
library(stringr)

#This sets the working directory
setwd("")

########################################This loads and preps the data 
#This is from 1.7
GROUPS <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This is from K10L2 3.4
MERGE_K10L2_RPM <- vroom::vroom("Tassel_TagTaxaDist_AllData_v7_v_K10L2.K10L2.RPM.txt")

#8 columns were classed as logicals because they have all NAs, these are the ones that had 0 coverage

#This is available in this repo under 4.3
BINS <- read_excel("Bins_NoOverlap_K10L2.table.xlsx")

#This drops any column that had an average coverage of 0, resulting in all values being NA
DT <- as.data.table(MERGE_K10L2_RPM)
DT <- DT[,which(unlist(lapply(DT, function(x)!all(is.na(x))))),with=F]
MERGE_K10L2_RPM_FILT_1 <- as.data.frame(DT)

#This saves the list of columns that were dropped
DROPPED <- setdiff(colnames(MERGE_K10L2_RPM), colnames(MERGE_K10L2_RPM_FILT_1))

#This removes any tag with a MAPQ less than 20
MERGE_K10L2_RPM_FILT <- subset(MERGE_K10L2_RPM_FILT_1, MAPQ >= 20)

###This selects only the controls
#This identifies all the index numbers for names marked as controls
col.num <- grep(".DC", colnames(MERGE_K10L2_RPM_FILT))
length(col.num)
#This appends the first 6 columns to that list
col.num <- append(col.num, values = c(1:6))
#This selects only the columns that appear in the above list
MERGE_K10L2_RPM_FILT_CONTROLS <- MERGE_K10L2_RPM_FILT[,sort(c(col.num))]

#This drops the Blanks
DROP <- grep("BLANK", colnames(MERGE_K10L2_RPM_FILT_CONTROLS) )
MERGE_K10L2_RPM_FILT_CONTROLS <- MERGE_K10L2_RPM_FILT_CONTROLS[,-DROP]
MISCLASS <-  c("W23_AB10-I.11.DC1", "W23_AB10-I.13.DC1", "W23_AB10-II.36.DC1", "W23_N10.14.DC1")
MERGE_K10L2_RPM_FILT_CONTROLS <- MERGE_K10L2_RPM_FILT_CONTROLS[,-c(which(colnames(MERGE_K10L2_RPM_FILT_CONTROLS) %in% MISCLASS))]

#This reduces the number of N10 samples included
GROUPS_CONTROLS_K10L2 <- subset(GROUPS, (GROUPS$Data_Source == "Dawe_Lab_1" | GROUPS$Data_Source == "Dawe_Lab_2") & (GROUPS$Ab10_Status == "K10L2"))

#This selects 5 N10 samples
GROUPS_CONTROLS_N10 <- subset(GROUPS, (GROUPS$Data_Source == "Dawe_Lab_1" | GROUPS$Data_Source == "Dawe_Lab_2") & (GROUPS$B_Chrom_Status == "No" & GROUPS$Ab10_Status == "N10"))

GROUPS_CONTROLS_N10 <- GROUPS_CONTROLS_N10[1:5,]

GROUPS_CONTROLS <- rbind(GROUPS_CONTROLS_K10L2, GROUPS_CONTROLS_N10)

#This identifies all the index numbers for names that appear in the GROUPS control only name fields
col.num <- which(colnames(MERGE_K10L2_RPM_FILT_CONTROLS) %in% GROUPS_CONTROLS$Name)
#This appends the first 6 columns to that list
col.num <- append(col.num, values = c(1:6))
#This selects only the columns that appear in the above list
MERGE_K10L2_RPM_FILT_CONTROLS <- MERGE_K10L2_RPM_FILT_CONTROLS[,sort(c(col.num))]


#this adds a bin value to each line in the MERGE_K10L2_RPM file 
MERGE_K10L2_RPM_FILT_CONTROLS$bin <- NA

i=1
x=1
for (i in 1:nrow(MERGE_K10L2_RPM_FILT_CONTROLS)) {
  print(i)
  for (x in 1:nrow(BINS)) {
    BIN_NUM <- BINS[x,4][[1]]
    START <- BINS[x,2][[1]]
    END <- BINS[x,3][[1]]
    ####CHANGE THIS NUMBER TO 3 WHEN YOU FIX THE PROBLEM WITH THE START
    if (MERGE_K10L2_RPM_FILT_CONTROLS[i,3] > START & MERGE_K10L2_RPM_FILT_CONTROLS[i,3] < END) {
      MERGE_K10L2_RPM_FILT_CONTROLS[i,ncol(MERGE_K10L2_RPM_FILT_CONTROLS)] <- BIN_NUM
    }
  }
}

#This outputs the bins as well as the tag coordinates to verify correct identification
test <- MERGE_K10L2_RPM_FILT_CONTROLS[,c(1:4,ncol(MERGE_K10L2_RPM_FILT_CONTROLS))]

#This function goes over each row and divides each value by the max in that row 
MinMax = function(xx) { sweep(xx, 1, apply(xx, 1, max), '/') }

#This sorts the data for plotting
SORT_temp1 <- subset(GROUPS, Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2")
SORT_temp2 <- SORT_temp1[order(SORT_temp1$Ab10_Status),]
SORT_temp3 <- SORT_temp2$Name
ADD <- c("Tag", "Chr", "Start", "End", "MAPQ", "Strand", "bin")
SORT_temp4 <- c(ADD, SORT_temp3)

#This removes the missing samples 
SORT <- SORT_temp4[! SORT_temp4 %in% setdiff(SORT_temp4, colnames(MERGE_K10L2_RPM_FILT_CONTROLS))]

#This manually alters SORT
SORT <- c(SORT[1:7], SORT[c(9,13)], SORT[c(8,10,11,12,14:25)])

#This orders the data frame
MERGE_K10L2_RPM_FILT_CONTROLS <- as.data.frame(MERGE_K10L2_RPM_FILT_CONTROLS[,SORT])

MERGE_K10L2_RPM_FILT_CONTROLS_TAGNUM <- MERGE_K10L2_RPM_FILT_CONTROLS

######## Tag Sums
#This sums the number of tag numbers across bins, but excludes the bin column from the summing
A_SUM = as.data.frame(apply(MERGE_K10L2_RPM_FILT_CONTROLS_TAGNUM[,c(8:ncol(MERGE_K10L2_RPM_FILT_CONTROLS_TAGNUM))], 2, function(xx) { by(xx, MERGE_K10L2_RPM_FILT_CONTROLS_TAGNUM$bin, sum)  }))

######## Tag Density
#This sums the number of tag numbers across bins, but excludes the bin column from the summing
A_DENS = as.data.frame(apply(MERGE_K10L2_RPM_FILT_CONTROLS_TAGNUM[,c(8:ncol(MERGE_K10L2_RPM_FILT_CONTROLS_TAGNUM))], 2, function(xx) { by(xx, MERGE_K10L2_RPM_FILT_CONTROLS_TAGNUM$bin, function(yy) { sum(yy > 0) }) }))

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

#Row 18 is an NA due to it always being 0. I am dropping it. The row lableed 18 is actually row 14
heat_data2 <- heat_data2[-c(14),]

#This generates a data frame with the Ab10 Status of all the lines after sorting them 
Name <- as.data.frame(colnames(heat_data2))
colnames(Name) <- "Name"
df_extra <- merge(Name, GROUPS)
df_col <- df_extra[,c("Name", "Ab10_Status")]
rownames(df_col) = df_col$Name

df_col_sub <- subset(df_col, Ab10_Status == "K10L2" | Ab10_Status == "N10")

#This orders the df_col rows to match the data columns
df_col_sub <- as.data.frame(df_col_sub[SORT[-c(1:7)],])

#This reorders the actual heat data to match above 
heat_data2 <- heat_data2[,SORT[-c(1:7)]]

#This finishes the df_col
df_col_sub <- as.data.frame(df_col_sub$Ab10_Status)
colnames(df_col_sub) <- c("Chr10 Type")
df_col_sub$`Chr10 Type` <- as.factor(df_col_sub$`Chr10 Type`)

#This line sets the rownames in the df as the column names in the heat map data. IT is very important
rownames(df_col_sub) = colnames(heat_data2)

#This creates a dataframe for the position data 
pos <- C_bed$start
df_row <- as.data.frame(pos)
df_row$feature <- " "
i=1
#This assigns Ab10 regions
for(i in 1:nrow(df_row)) {
  VALUE <- df_row[i,1]
  df_row[i,2] <- ifelse(VALUE >= 0 & VALUE <= 1381015, 'Shared Region', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 1383651 & VALUE <= 7920490, 'TR-1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 10652268 & VALUE <= 19464951, 'TR-1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 19450000 & VALUE <= 25225000, 'Shared Region', df_row[i,2])
}
#Trkin is too small for this loop to assign it well for some reason. It should be number 11
df_row[11,2] <- "Trkin"


df_row <- as.data.frame(df_row$feature)
rownames(df_row) = rownames(C_bed)
colnames(df_row) <- "Feature"
df_row$Feature <- as.factor(df_row$Feature)


heat_data2 <- MinMax(heat_data2)

colors <- list(`Chr10 Type` = c("N10" ="grey70", "K10L2" = "black"), Feature = c('TR-1' = "#56B4E9", 'Trkin' = "#0072B2", 'Shared Region' = "#E69F00", " " = "white" ))

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_K10L2_recolor_publish.pdf", height = 5, width = 3)
pheatmap(heat_data2, scale = "none", show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols =FALSE, clustering_method = "ward.D", treeheight_row = 0, treeheight_col = 0, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), border_color = NA, annotation_col = df_col_sub, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_K10L2_recolor.pdf", height = 5, width = 5)
pheatmap(heat_data2, scale = "none", show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE, treeheight_row = 0, treeheight_col = 0, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), border_color = NA, annotation_col = df_col_sub, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_K10L2_ColNames_recolor.pdf", height = 5, width = 5)
pheatmap(heat_data2, scale = "none", show_rownames = FALSE, show_colnames = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, treeheight_row = 0, treeheight_col = 0, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), border_color = NA, annotation_col = df_col_sub, fontsize_col = 2, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_K10L2_ColClust_recolor.pdf", height = 5, width = 5)
pheatmap(heat_data2, scale = "none", show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = TRUE, clustering_method = "ward.D", treeheight_row = 0, treeheight_col = 50, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500),  border_color = NA, annotation_col = df_col_sub, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_K10L2_ColClust_ColNames_recolor.pdf", height = 5, width = 10)
pheatmap(heat_data2, scale = "none", show_rownames = FALSE, show_colnames = TRUE, cluster_rows = FALSE, cluster_cols = TRUE, clustering_method = "ward.D", treeheight_row = 0, treeheight_col = 50, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500),  border_color = NA, annotation_col = df_col_sub, fontsize_col = 3, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_K10L2_ColClust_RowClust_recolor.pdf", height = 5, width = 10)
pheatmap(heat_data2, scale = "none", show_rownames = FALSE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D", treeheight_row = 50, treeheight_col = 50, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), border_color = NA, annotation_col = df_col_sub, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_B73Ab10v2_TagIndexMinMax_K10L2_ColClust_RowClust_RowNames_recolor.pdf", height = 5, width = 10)
pheatmap(heat_data2, scale = "none", show_rownames = TRUE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D", treeheight_row = 50, treeheight_col = 50, annotation_row = df_row, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), border_color = NA, annotation_col = df_col_sub, fontsize_row = 3, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()
