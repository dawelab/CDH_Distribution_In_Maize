#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)
IT <-args[1]
print(IT)

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

#This sets the working directory where all files being loaded in are located
setwd("")

#This loads the original groups file
#from 1.7
GROUPS <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This loads in the groups file with modified names
#from 6.1.1
DF <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table")

#This loads in the BINS file
#from 4.3
BINS <- read_excel("Bins_NoOverlap_K10L2.table.xlsx")

#From 6.2.3
POS <- vroom::vroom("Kmeans_Exp_Groups_MixedControls_All.csv")

#This function goes over each row and divides each value by the max in that row 
MinMax = function(xx) { sweep(xx, 1, apply(xx, 1, max), '/') }

#This loads in the combination of the tag sums and tag density
#This is from 6.1.4
B_COMB <- vroom::vroom("K10L2_TagSumPlusTagDensAllSamples.csv")

#This scales the data so that the min is always 0 and the max is always 1. 
C <- MinMax(B_COMB)

#These lines add bins values back in, they were dropped in the ASum and ADens functions. 
C$bin <- rownames(C)
C_bed <- merge(BINS, C, by="bin")

#This drops the identifier columns to create a numeric matrix, but saves the bin value as the row number
rownames(C_bed) <- C_bed$bin
heat_data <- C_bed[,-c(1:4)]

#This selects only unique entries, columns 2 and 3 include the extra Romero Navarro names to differentiate each replicate and need to be removed to find unique lines
DF_EXTRA <- unique(DF[,-c(2:3)])

#This selects only control data, but excludes K10L2 samples and B73 samples which need to be processed separatly. It also drops control samples known to be misclassified. 
DF_EXTRA_2 <- subset(DF_EXTRA, (Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2") & Ab10_Status != "Ab10_Unknown" & Ab10_Status != "Ab10-II" & Name != "W23_AB10-I.11.DC1"& Name != "W23_AB10-I.13.DC1" & Name != "W23_AB10-II.36.DC1" & Name != "W23_AB10-II.37.DC1" & Name != "W23_N10.14.DC1" & Name != "PI-483314_K10L2.1.DC2")

#This isolates only the controls 
heat_data_control <- heat_data[,DF_EXTRA_2$Name]

#This isolates only the experimental samples
DF_EXTRA_3 <- subset(DF_EXTRA, Data_Source != "Dawe_Lab_1" & Data_Source != "Dawe_Lab_2")

#This selects only samples that were called N10 is the Ab10 model
#This selects Ab10 positive experimental samples
POS <- POS[,-c(127)]

POS$All <- rowSums(POS[,-c(1)])
N10_sub <- subset(POS, All < 125)
N10_POS <- N10_sub$Name

#This selects only N10 experimental samples
DF_EXTRA_4 <- DF_EXTRA_3[which(DF_EXTRA_3$Name %in% N10_POS),]

#These lines generate a dataframe that will be used to annotate columns by chromosome 10 haplotype and data source of origin. 
df_col <- DF_EXTRA[,c("Name", "Ab10_Status", "Data_Source")]
colnames(df_col) <- c("Name", "Chr10 Type", "Data Source")
df_col$`Chr10 Type` <- as.factor(df_col$`Chr10 Type`)
df_col$`Data Source` <- as.factor(df_col$`Data Source`)

#This creates a dataframe that will be used to annotate the rows by Ab10 feature
pos <- C_bed$start
df_row <- as.data.frame(pos)
df_row$feature <- " "

#This assigns Ab10 regions
for(i in 1:nrow(df_row)) {
  VALUE <- df_row[i,1]
  df_row[i,2] <- ifelse(VALUE >= 1383651 & VALUE <= 7920490, 'TR1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 10652268 & VALUE <= 19464951, 'TR1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 10243298 & VALUE <= 10319813, 'trkin', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 19450000 & VALUE <= 25225000, 'Shared region', df_row[i,2])
}
#Trkin is too small for this loop to assign it well for some reason. It should be number 11
df_row[11,2] <- "trkin"

df_row <- as.data.frame(df_row$feature)
rownames(df_row) = rownames(C_bed)
colnames(df_row) <- "Feature"
df_row$Feature <- as.factor(df_row$Feature)

#This determines the column annotation colors 
columntype_colors <- c("K10L2" = "brown","Ab10-I" = "black","Ab10-II" = "grey20","Ab10-III" =  "grey40", "N10" ="grey80", "Unknown" = "white")

columndatasource_colors <- c("Dawe_Lab_1" = "gold", "Dawe_Lab_2" = "limegreen", "Romero-Navarro_etal_2017" = "#008080", "Swarts_etal_2017" = "#310062", "Romay_etal_2013" = "mediumorchid", "NA" = "white")

row_colors <- c('TR1' = "#56B4E9", 'trkin' = "#0072B2", 'Shared region' = "#E69F00")

#This determines the grid colors
grid_color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500)

#This creates a data frame to add the Kmeans groupings too
KMeans_Groups_All <- data.frame(Name = c(NA), Group = c(NA), Ab10_Status = c(NA), Correct = c(NA))
KMeans_Groups_All <- KMeans_Groups_All[-c(1),]

#This creates the final data frame
Final <- data.frame(Name = colnames(heat_data))

#This writes out a data drame to ensure all the controls are being classed correctly 
QC_Controls <- data.frame(Iteration=c(NA), Round=c(NA), Subround=c(NA), Controls_Correct=c(NA))
QC_Controls <- QC_Controls[-c(1),]

###############################################################################

#heat_data_control is the matrix containing the heat data
#KMeans_Groups_All is the temporary data frame to add each run of all controls to
#j is the iteration
for(j in 1:5) {
  #This is the iteration
  print(paste("This is iteration number", j, sep=" "))
  #This creates a vector with all of the experimental sample names
  TEMP_EXP <- DF_EXTRA_4$Name
  #These lines determine the number of iterations needed to cluster each experimental sample, while only introduces 25% experimental samples with all others being control samples. q is the rounded down whole number, r is the remainder. s_EXP is a vector repeating 25% of the number of controls q times and then adding r to 25% the number of controls for the last number. 
  q <- floor((length(TEMP_EXP))/(round((ncol(heat_data_control)*.25))))
  
  r <- length(TEMP_EXP)-(q*(round((ncol(heat_data_control)*.25))))
  
  s_EXP <- c(rep(round((ncol(heat_data_control)*.25)),q-1), (round((ncol(heat_data_control)*.25))+r))
  
  #This creates a data frame to add the Kmeans groupings too
  KMeans_Groups_All <- data.frame(Name = c(NA), Group = c(NA), Ab10_Status = c(NA), K10L2_Status = c(NA), Correct = c(NA))
  KMeans_Groups_All <- KMeans_Groups_All[-c(1),]
  
  #This randomly selects the appropriate number of experimental samples as determined above, adds them to the controls, and writes out the file. It then removes the selected experimental samples from the pool and repeat's the process. In this way each experimental sample is only selected once. 
  
  for(i in 1:q) {
    #This pulls the appropriate number of experimental samples for this iteration
    n_EXP <- s_EXP[i]
    #This randomly selects n_EXP experimental samples from the vector or all experimental samples without replacement by index value. 
    samp_EXP <- sample(1:length(TEMP_EXP), n_EXP, replace = FALSE)
    #This pulls the randomly selected index values from the TEMP_EXP vector, returning actual sample names and subsets the large heat_data matrix to those columns by name. 
    heat_data_exp <- heat_data[,c(TEMP_EXP[samp_EXP])]
    #This adds the experimental samples to the control samples
    heat_data_samp <- cbind(heat_data_control, heat_data_exp)
    #This writes out the file
    write.table(heat_data_samp, file = paste("K10L2_Model/Experimental_K10L2N10Model_v2/SubSamples/K10L2N10_Exp_heat_data_samp_MixedControls", IT, i, "txt", sep ="."), row.names = FALSE, quote = FALSE, sep = "\t")
    #This removes the experimental samples from the list of available experimental samples so that they are not chosen again. 
    TEMP_EXP <- TEMP_EXP[-samp_EXP]
  }

  
  k <- 1
  while(k <= q) {
    print(k)
    heat_data_samp <- read.table(paste("K10L2_Model/Experimental_K10L2N10Model_v2/SubSamples/K10L2N10_Exp_heat_data_samp_MixedControls", IT, k, "txt", sep ="."), header = TRUE,  sep = "\t", check.names=FALSE)
    df_col_samp <- df_col[which(df_col$Name %in% colnames(heat_data_samp)),]
    heat_data_samp <- heat_data_samp[,df_col_samp$Name]
    #This groups by Kmeans
    group = kmeans(t(heat_data_samp[c(7:10),]), centers = 2)$cluster
    
    #This interprets and writes out the kmeans 
    KMeans_Groups <-  data.frame(Name=colnames(heat_data_samp), Group=group)
    
    KMeans_Groups_2 <- merge(KMeans_Groups, GROUPS, by = "Name")
    KMeans_Groups_3 <- KMeans_Groups_2[,c("Name", "Group.x", "Ab10_Status")]
    colnames(KMeans_Groups_3) <- c("Name", "Group", "Ab10_Status")
    
    KMeans_Groups_3$K10L2_Status <- KMeans_Groups_3$Ab10_Status
    
    KMeans_Groups_3$K10L2_Status <-  gsub("Ab10-III", "K10L2", KMeans_Groups_3$K10L2_Status)
    KMeans_Groups_3$K10L2_Status <-  gsub("Ab10-I", "K10L2", KMeans_Groups_3$K10L2_Status)
    
    #This calculates the percent of each Kmeans group that is called Ab10 or N10
    Percent <- KMeans_Groups_3 %>%
      group_by(Group) %>% 
      count(K10L2_Status) %>%
      group_by(Group) %>%
      mutate(percent = n/sum(n))
    
    #This determines the predominant call for each Kmeans groups
    #This does group 1
    if(length(grep(1, Percent$Group)) > 1) {
      sub <- subset(Percent, Percent$Group == 1)
      sub2 <- subset(sub, K10L2_Status != "Unknown")
      if(nrow(sub2) == 1) {
        KMeans_Groups_3$Group <- gsub("\\b1\\b", sub[1,2], KMeans_Groups_3$Group, perl = TRUE)
      } else {
        KMeans_Groups_3$Group <- gsub("\\b1\\b", paste("Mixed", nrow(sub2), sub2[1,2], round(sub2[1,4], digits=2), sub2[2,2], round(sub2[2,4], digits=2), sep="_"), KMeans_Groups_3$Group, perl = TRUE)
      }
    } else{
      sub <- subset(Percent, Percent$Group == 1)
      KMeans_Groups_3$Group <- gsub("\\b1\\b", sub[1,2], KMeans_Groups_3$Group, perl = TRUE)
    }
    
    #This does group 2
    if(length(grep(2, Percent$Group)) > 1) {
      sub <- subset(Percent, Percent$Group == 2)
      sub2 <- subset(sub, K10L2_Status != "Unknown")
      if(nrow(sub2) == 1) {
        KMeans_Groups_3$Group <- gsub("\\b2\\b", sub[1,2], KMeans_Groups_3$Group, perl = TRUE)
      } else {
        KMeans_Groups_3$Group <- gsub("\\b2\\b", paste("Mixed", nrow(sub2), sub2[1,2], round(sub2[1,4], digits=2), sub2[2,2], round(sub2[2,4], digits=2), sep="_"), KMeans_Groups_3$Group, perl = TRUE)
      }
    } else{
      sub <- subset(Percent, Percent$Group == 2)
      KMeans_Groups_3$Group <- gsub("\\b2\\b", sub[1,2], KMeans_Groups_3$Group, perl = TRUE)
    }
    
    #This determines control correctness 
    KMeans_Groups_4 <- subset(KMeans_Groups_3, Ab10_Status != "Unknown")
    KMeans_Groups_4$Correct <- ifelse(KMeans_Groups_4$Group == KMeans_Groups_4$K10L2_Status, 1, 0)
    
    #This loop checks that all controls were called correctly. If they are it proceeds with the loop, if not it prints a message and redoes that same number again
    if(sum(KMeans_Groups_4$Correct) == length(KMeans_Groups_4$Correct)) {
      
      #This reports the correct controls
      line <- data.frame(Iteration=c(IT), Round=c(j), Subround=c(k), Controls_Correct=sum(KMeans_Groups_4$Correct)/nrow(KMeans_Groups_4))
      QC_Controls <<- rbind(QC_Controls, line)
      
      #This stores the calls for this run 
      KMeans_Groups_5 <- subset(KMeans_Groups_3, Ab10_Status == "Unknown")
      KMeans_Groups_All <<- rbind(KMeans_Groups_All, KMeans_Groups_5)
      
      #This makes the heatmap annotations
      ha_col = HeatmapAnnotation(`Chr10 Type` = df_col_samp$`Chr10 Type`,
                                 `Data Source` = df_col_samp$`Data Source`,
                                 col = list(`Chr10 Type` = columntype_colors, `Data Source` = columndatasource_colors))
      
      ha_row = rowAnnotation(Feature = df_row$Feature,
                             col = list(Feature = row_colors),
                             na_col = "white")
      
      heat_data_samp_rescale <- MinMax(heat_data_samp)
      
      pdf(file = paste("K10L2_Model/Experimental_K10L2N10Model_v2/Images/K10L2N10_HeatMap_MixedControls", IT, j, k,"pdf", sep="."))
      a <- Heatmap(as.matrix(heat_data_samp_rescale),
                   #plotting
                   col = grid_color,
                   column_title = "Tag Index Across K10L2 Haplotype",
                   column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                   row_title = "K10L2 Haplotype",
                   border_gp = gpar(col = "black", lty = 1),
                   
                   #clustering
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   column_split = group,
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
      
      #This adds one to the value of k
      k <<- k+1
    } else {
      print(paste("This iteration had incorrectly grouped controls and was redone", IT, ".", j, ".", k))
    }
  }
  
  KMeans_Groups_All_v2 <- KMeans_Groups_All
  colnames(KMeans_Groups_All_v2) <- c("Name", paste("Group", IT, j, sep="_"), "Ab10_Status")
  #This adds this rounds data to the output dataframe
  KMeans_Groups_All_v2_sub <- KMeans_Groups_All_v2[,c(1,2)]
  Final <<- merge(Final,KMeans_Groups_All_v2_sub, by="Name")
}

fwrite(Final, file = paste("K10L2_Model/Experimental_K10L2N10Model_v2/KmeansGroups_Files/K10L2_Kmeans_Exp_Groups_MixedControls", IT, "table", sep="."))

fwrite(QC_Controls, file = paste("K10L2_Model/Experimental_K10L2N10Model_v2/QualityControl/K10L2_QualityControl_KMeansGroupingControls_MixedControls", IT, "table", sep="."))
