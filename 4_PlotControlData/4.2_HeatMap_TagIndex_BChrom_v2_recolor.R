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

MERGE_BChrom_RPM<- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.BChrom.RPM.txt")

#7 columns were classed as logicals because they have all NAs, these are the ones that had 0 coverage

BINS <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Bins_NoOverlap_BChrom.table.xlsx")

#This drops any column that had an average coverage of 0, resulting in all values being NA
DT <- as.data.table(MERGE_BChrom_RPM)
DT <- DT[,which(unlist(lapply(DT, function(x)!all(is.na(x))))),with=F]
MERGE_BChrom_RPM<- as.data.frame(DT)

#This removes any tag with a MAPQ less than 20
MERGE_BChrom_RPM_FILT <- subset(MERGE_BChrom_RPM, MAPQ >= 20)

###This selects only the controls
#This selects only the control data
GROUPS_CONTROLS <- subset(GROUPS, GROUPS$Data_Source == "Dawe_Lab_1" | GROUPS$Data_Source == "Dawe_Lab_2")

#This identifies all the index numbers for names that appear in the GROUPS control only name fields
col.num <- which(colnames(MERGE_BChrom_RPM_FILT) %in% GROUPS_CONTROLS$Name)
#This appends the first 6 columns to that list
col.num <- append(col.num, values = c(1:6))
#This selects only the columns that appear in the above list
MERGE_BChrom_RPM_FILT_CONTROLS <- MERGE_BChrom_RPM_FILT[,sort(c(col.num))]

#This drops lines that are are not relevant to identifying B Chromosomes
DROP_1 <- c("B73_N10.1.DC2", "B73_N10.2.DC2", "B73_N10.3.DC2", "BLANK.1.DC1", "BLANK.2.DC1", "BLANK.1.DC2", "NSL-2833_B-Chrom.2.DC2", "B542C_L289_B-Chrom.1.DC2")
K10L2 <- subset(GROUPS, Ab10_Status == "K10L2")
DROP_2 <- c(DROP_1, K10L2$Name)
DROP_3 <- which(colnames(MERGE_BChrom_RPM_FILT_CONTROLS)  %in% DROP_2)
MERGE_BChrom_RPM_FILT_CONTROLS <- MERGE_BChrom_RPM_FILT_CONTROLS[,-DROP_3]

colnames(MERGE_BChrom_RPM_FILT_CONTROLS)

#This reduces the number of N10 samples included
GROUPS_CONTROLS_B <- subset(GROUPS, (GROUPS$Data_Source == "Dawe_Lab_1" | GROUPS$Data_Source == "Dawe_Lab_2") & (GROUPS$B_Chrom_Status == "Yes"))

#This selects 5 N10 samples
GROUPS_CONTROLS_N10 <- subset(GROUPS, (GROUPS$Data_Source == "Dawe_Lab_1" | GROUPS$Data_Source == "Dawe_Lab_2") & (GROUPS$B_Chrom_Status == "No" & GROUPS$Ab10_Status == "N10"))

GROUPS_CONTROLS_N10 <- GROUPS_CONTROLS_N10[1:5,]

GROUPS_CONTROLS <- rbind(GROUPS_CONTROLS_B, GROUPS_CONTROLS_N10)

#This identifies all the index numbers for names that appear in the GROUPS control only name fields
col.num <- which(colnames(MERGE_BChrom_RPM_FILT_CONTROLS) %in% GROUPS_CONTROLS$Name)
#This appends the first 6 columns to that list
col.num <- append(col.num, values = c(1:6))
#This selects only the columns that appear in the above list
MERGE_BChrom_RPM_FILT_CONTROLS <- MERGE_BChrom_RPM_FILT_CONTROLS[,sort(c(col.num))]

#this adds a bin value to each line in the MERGE_BChrom_RPMfile 
MERGE_BChrom_RPM_FILT_CONTROLS$bin <- NA
i=1
x=1
for (i in 1:nrow(MERGE_BChrom_RPM_FILT_CONTROLS)) {
  print(i)
  for (x in 1:nrow(BINS)) {
    BIN_NUM <- BINS[x,4][[1]]
    START <- BINS[x,2][[1]]
    END <- BINS[x,3][[1]]
    if (MERGE_BChrom_RPM_FILT_CONTROLS[i,3] > START & MERGE_BChrom_RPM_FILT_CONTROLS[i,3] < END) {
      MERGE_BChrom_RPM_FILT_CONTROLS[i,ncol(MERGE_BChrom_RPM_FILT_CONTROLS)] <- BIN_NUM
    }
  }
}

#This outputs the bins as well as the tag coordinates to verify correct identification
test <- MERGE_BChrom_RPM_FILT_CONTROLS[,c(1:4,ncol(MERGE_BChrom_RPM_FILT_CONTROLS))]

#This function goes over each row and divides each value by the max in that row 
MinMax = function(xx) { sweep(xx, 1, apply(xx, 1, max), '/') }

#This is being a problem and I'm not sure why 
SORT_temp1 <- subset(GROUPS, Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2")
SORT_temp2 <- SORT_temp1[order(SORT_temp1$B_Chrom_Status),]
SORT_temp3 <- SORT_temp2$Name
ADD <- c("Tag", "Chr", "Start", "End", "MAPQ", "Strand", "bin")
SORT_temp4 <- c(ADD, SORT_temp3)

#This removes the missing samples 
SORT <- SORT_temp4[! SORT_temp4 %in% setdiff(SORT_temp4, colnames(MERGE_BChrom_RPM_FILT_CONTROLS))]
#By default the B yes' are on the right I want them to be on the left in the final plot
SORT <- c(SORT[1:7], rev(SORT[-c(1:7)]))

#This orders the data frame
MERGE_BChrom_RPM_FILT_CONTROLS <- as.data.frame(MERGE_BChrom_RPM_FILT_CONTROLS[,SORT])

######## Tag Sums
#This sums the number of tag numbers across bins, but excludes the bin column from the summing
A_SUM = as.data.frame(apply(MERGE_BChrom_RPM_FILT_CONTROLS[,c(8:ncol(MERGE_BChrom_RPM_FILT_CONTROLS))], 2, function(xx) { by(xx, MERGE_BChrom_RPM_FILT_CONTROLS$bin, sum)  }))

######## Tag Density
#This sums the number of tag numbers across bins, but excludes the bin column from the summing
A_DENS = as.data.frame(apply(MERGE_BChrom_RPM_FILT_CONTROLS[,c(8:ncol(MERGE_BChrom_RPM_FILT_CONTROLS))], 2, function(xx) { by(xx, MERGE_BChrom_RPM_FILT_CONTROLS$bin, function(yy) { sum(yy > 0) }) }))

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
heat_data <- as.matrix(C_bed[,-c(1:4)])

#This reorders the columns 
heat_data <- heat_data[,c(1,3,2,4,6,10,11,12,13,5,7,8,9,14,15,16,17,18)]

#write.table(cols, file = "column_names_dens_RPM.csv", row.names = FALSE)
#This generates a data frame with the Ab10 Status of all the lines after sorting them 
Name <- as.data.frame(colnames(heat_data))
colnames(Name) <- "Name"
df_extra <- merge(Name, GROUPS)
df_col <- df_extra[,c("Name", "B_Chrom_Status")]
rownames(df_col) = df_col$Name

#This sorts the data frame
df_col <- df_col[SORT[-c(1:7)],]

rownames(df_col) <- df_col$Name
df_col <- as.data.frame(df_col$B_Chrom_Status)
colnames(df_col) <- c("B Chromosome")
df_col$`B Chromosome` <- as.factor(df_col$`B Chromosome`)

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

#df_row[1,2] <- 'B Short Arm'

df_row <- as.data.frame(df_row$feature)
rownames(df_row) = rownames(C_bed)
colnames(df_row) <- "Feature"
df_row$Feature <- ifelse(is.na(df_row$Feature)==TRUE, " ",df_row$Feature )
df_row$Feature <- factor(df_row$Feature, levels= c('B Repeat', 'knob180', 'Distal Heterochromatin', ' '))

rownames(df_col) <- colnames(heat_data)

colors <- list(`B Chromosome` = c("Yes" = "black","No" = "grey70"), Feature = c('B Repeat' = "lightpink2", 'knob180' = "#D55E00",'Distal Heterochromatin' =  "#009e73", " " = "white"))

pdf("HeatMap_TagSumsSqrtTagDens_MinMax_ColClust_BChrom_recolor_publish.pdf", height = 5, width = 3.5)
pheatmap(heat_data, scale = "none", show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE, treeheight_row = 0, treeheight_col = 0, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), annotation_row = df_row, annotation_col = df_col, border_color = NA, main="B Chromosome\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_TagSumsSqrtTagDens_MinMax_ColClust_BChrom_recolor.pdf", height = 5, width = 5)
pheatmap(heat_data, scale = "none", show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = TRUE, treeheight_row = 0, treeheight_col = 50, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), annotation_row = df_row, annotation_col = df_col, main="B Chromosome\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_TagSumsSqrtTagDens_MinMax_ColClust_ColNames_BChrom_recolor.pdf", height = 5, width = 10)
pheatmap(heat_data, scale = "none", show_rownames = FALSE, show_colnames = TRUE, cluster_rows = FALSE, cluster_cols = TRUE, treeheight_row = 0, treeheight_col = 50, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), annotation_row = df_row, annotation_col = df_col, fontsize_col = 3, main="B Chromosome\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_TagSumsSqrtTagDens_MinMax_ColClust_RowClust_BChrom_recolor.pdf", height = 5, width = 10)
pheatmap(heat_data, scale = "none", show_rownames = FALSE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, treeheight_row = 50, treeheight_col = 50, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), annotation_row = df_row, annotation_col = df_col, main="B Chromosome\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_TagSumsSqrtTagDens_MinMax_ColClust_RowClust_RowNames_BChrom_recolor.pdf", height = 5, width = 10)
pheatmap(heat_data, scale = "none", show_rownames = TRUE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, treeheight_row = 50, treeheight_col = 50, annotation_colors= colors, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500), annotation_row = df_row, annotation_col = df_col, fontsize_row = 3, main="B Chromosome\nTag Index in 1MB bins")
dev.off()

