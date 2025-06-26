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
library(ComplexHeatmap)
library(MASS)
install.packages("randomForest")
library(randomForest)
install.packages("caret")
library(caret)

########################################This loads and preps the data 
#This file is from 1.7
GROUPS <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This file is from 3.4
MERGE_Ab10Hap_RPM_FILT_2_FIX <- vroom::vroom("BWAaln_All_v_B73-Ab10_BChrom.Ab10.RPM.RNMean.table")

#This is from 6.1.1
DF <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table")

#This loads in the BINS file
#This is from 4.1
BINS <- read_excel("Bins_NoOverlap.table.xlsx")

#This is from 6.2.3 and 4
POS <- vroom::vroom("Kmeans_Exp_Groups_MixedControls_All.csv")

#This function goes over each row and divides each value by the max in that row 
MinMax = function(xx) { sweep(xx, 1, apply(xx, 1, max), '/') }

######## Tag Sums
#This sums the number of tag numbers across bins, but excludes the bin column from the summing

A_SUM = as.data.frame(apply(MERGE_Ab10Hap_RPM_FILT_2_FIX[,-c(1:7)], 2, function(xx) { by(xx, MERGE_Ab10Hap_RPM_FILT_2_FIX$bin, sum)  }))

fwrite(A_SUM, file ="Ab10_TagSums_AllSamples.csv")

######## Tag Density
#This sums the number of tag numbers across bins, but excludes the bin column from the summing
A_DENS = as.data.frame(apply(MERGE_Ab10Hap_RPM_FILT_2_FIX[,-c(1:7)], 2, function(xx) { by(xx, MERGE_Ab10Hap_RPM_FILT_2_FIX$bin, function(yy) { sum(yy > 0) }) }))

fwrite(A_SUM, file ="Ab10_TagDens_AllSamples.csv")

#This takes the square root of each value

B_DENS <- sqrt(A_DENS)

#This combines the density and sum data
B_COMB <- A_SUM+B_DENS

fwrite(B_COMB, file ="Ab10_TagSumPlusTagDensAllSamples.csv")

#This scales the data so that the min is always 0 and the max is always 1. 
C <- MinMax(B_COMB)

#This adds bins back
C$bin <- rownames(C)

C_bed <- merge(BINS, C, by="bin")

#This drops the identifier columns to create a numeric matrix
rownames(C_bed) <- C_bed$bin
heat_data <- C_bed[,-c(1:4)]


#This selects only unique entries, columns 2 and 3 include the extra Romero Navarro names to differentiate each replicate and need to be removed to find unique lines
DF_EXTRA <- unique(DF[,-c(2:3)])
DF_EXTRA_2 <- subset(DF_EXTRA, (Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2") & Ab10_Status != "K10L2" & Ab10_Status != "Ab10_Unknown" & Name != "B73_N10.1.DC2" & Name != "B73_N10.2.DC2" & Name != "B73_N10.3.DC2" & Name != "W23_AB10-I.11.DC1"& Name != "W23_AB10-I.13.DC1" & Name != "W23_AB10-II.36.DC1" & Name != "W23_N10.14.DC1" & Name != "Ab10-I-Jal_het.1.DC2" & Name != "Ab10-I-Jal_het.2.DC2" & Name != "Ab10-I-Jal_het.3.DC2" & Name != "Ab10-I-Jal_het.4.DC2" & Name != "Ab10-I-Jal_het.5.DC2")

#This subsets only the Ab10 controls 
heat_data_control <- heat_data[,DF_EXTRA_2$Name]

#This selects Ab10 positive experimental samples
POS <- POS[,-c(127)]
POS$All <- rowSums(POS[,-c(1)])
POS_sub <- subset(POS, All < 125 & All > 0)
Ab10_AMBIG <- POS_sub$Name

#This subsets the heat data to the high confidence type calls 
heat_data_Ambig <- heat_data[, which(colnames(heat_data) %in% Ab10_AMBIG)]

#This merges the controls and the Ambig
heat_data_controlAmbig <- cbind(heat_data_control, heat_data_Ambig)

#Select only the relevant columns
df_col <- DF_EXTRA[,c("Name", "Ab10_Status", "Data_Source")]
colnames(df_col) <- c("Name", "Chr10 Type", "Data Source")
df_col$`Chr10 Type` <- as.factor(df_col$`Chr10 Type`)
df_col$`Data Source` <- as.factor(df_col$`Data Source`)

#This creates the df_col dataframe 
df_col_controlAmbig <- df_col[which(df_col$Name %in% colnames(heat_data_controlAmbig)),]

#This section reorders the data to match the heat data
df_col_controlAmbig <- df_col_controlAmbig[order(df_col_controlAmbig$`Chr10 Type`),]
heat_data_controlAmbig <- heat_data_controlAmbig[,c(df_col_controlAmbig$Name)]

#This creates a dataframe for the position data 
pos <- C_bed$start
df_row <- as.data.frame(pos)
df_row$feature <- " "

#This assigns Ab10 regions
for(i in 1:nrow(df_row)) {
  VALUE <- df_row[i,1]
  df_row[i,2] <- ifelse(VALUE >= 142472000 & VALUE <= 146699300, 'TR1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 148964528 & VALUE <= 149082763, 'trkin', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 150656000 & VALUE <= 153145000, 'TR1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 152050000 & VALUE <= 156350000, 'Shared region', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 158250000 & VALUE <= 166820000, 'Shared region', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 157485200 & VALUE <= 159356550, 'TR1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 174433450 & VALUE <= 182846100, 'knob 180', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 189326066 & VALUE <= 190330226, 'kindr', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 190307761 & VALUE <= 191107857, 'kin10-like', df_row[i,2])
}

df_row <- as.data.frame(df_row$feature)
rownames(df_row) = rownames(C_bed)
colnames(df_row) <- "Feature"
df_row$Feature <- as.factor(df_row$Feature)

#This determines the column annotation colors 
#This determines the grid colors
grid_color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500)

columntypetrue_colors <- c("Ab10-I" = "black","Ab10-II" = "grey20","Ab10-III" =  "grey40", "N10" ="grey80", "Unknown" = "white", "N10" = "beige")

columndatasource_colors <- c("Dawe_Lab_1" = "gold", "Dawe_Lab_2" = "limegreen", "Romero-Navarro_etal_2017" = "#008080", "Swarts_etal_2017" = "#310062", "Romay_etal_2013" = "mediumorchid", "NA" = "white")

rowfeature_colors <- c('TR1' = "#56B4E9", 'trkin' = "#0072B2", 'Shared region' = "#E69F00", 'knob 180' = "#D55E00", 'kindr' = "#CC79A7", 'kin10-like' =  "#009E73")


#This makes the heatmap annotations
ha_col = HeatmapAnnotation(`Chr10 Type` = df_col_controlAmbig$`Chr10 Type`,
                           `Data Source` = df_col_controlAmbig$`Data Source`,
                           col = list(`Chr10 Type` = columntypetrue_colors, `Data Source` = columndatasource_colors))

ha_row = rowAnnotation(Feature = df_row$Feature,
                       col = list(Feature = rowfeature_colors),
                       na_col = "white")

pdf(file = paste("KMeansAb10N10Ambig_HeatMap", j, "pdf", sep="."), height=10, width=10)

a <- Heatmap(as.matrix(heat_data_controlAmbig),
             #plotting
             col = grid_color,
             column_title = "Tag Index Across Ab10 Haplotype",
             column_title_gp = gpar(fontsize = 20, fontface = "bold"),
             row_title = "Ab10 Haplotype", 
             border_gp = gpar(col = "black", lty = 1),
             
             #clustering
             cluster_rows = FALSE, 
             cluster_columns = FALSE, 
             column_split = df_col_controlAmbig$`Chr10 Type`,
             column_gap = unit(5, "mm"),
             
             #annotations
             name="Tag Index",
             top_annotation = ha_col,
             left_annotation = ha_row, 
             show_row_names= TRUE,
             show_column_names = TRUE, 
             row_names_gp = grid::gpar(fontsize = 5),
             column_names_gp = grid::gpar(fontsize = 0.5)
)
print(a)
dev.off()

