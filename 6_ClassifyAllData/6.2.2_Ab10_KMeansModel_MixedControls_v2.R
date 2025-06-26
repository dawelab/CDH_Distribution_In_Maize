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
#This is from 1.7
GROUPS <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This loads in the groups file with modified names
#This is from 6.1.1
DF <- vroom::vroom("Ab10_Model/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table")

#This loads in the BINS file
#This is from 4.1
BINS <- read_excel("Bins_NoOverlap.table.xlsx")

#This function goes over each row and divides each value by the max in that row 
MinMax = function(xx) { sweep(xx, 1, apply(xx, 1, max), '/') }

#This loads in the combination of the tag sums and tag density
B_COMB <- vroom::vroom("Ab10_Model/Experimental_Ab10N10Model/Ab10_TagSumPlusTagDensAllSamples.csv")

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
DF_EXTRA_2 <- subset(DF_EXTRA, (Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2") & Ab10_Status != "K10L2" & Ab10_Status != "Ab10_Unknown" & Name != "B73_N10.1.DC2" & Name != "B73_N10.2.DC2" & Name != "B73_N10.3.DC2" & Name != "W23_AB10-I.11.DC1"& Name != "W23_AB10-I.13.DC1" & Name != "W23_AB10-II.36.DC1" & Name != "W23_N10.14.DC1")

#This isolates only the experimental samples
DF_EXTRA_3 <- subset(DF_EXTRA, Data_Source != "Dawe_Lab_1" & Data_Source != "Dawe_Lab_2")

#This isolates only the controls 
heat_data_control <- heat_data[,DF_EXTRA_2$Name]

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
columntype_colors <- c("Ab10-I" = "black","Ab10-II" = "grey20","Ab10-III" =  "grey40", "N10" ="grey80", "Unknown" = "white")

columndatasource_colors <- c("Dawe_Lab_1" = "gold", "Dawe_Lab_2" = "limegreen", "Romero-Navarro_etal_2017" = "#008080", "Swarts_etal_2017" = "#310062", "Romay_etal_2013" = "mediumorchid", "NA" = "white")

row_colors <- c('TR1' = "#56B4E9", 'trkin' = "#0072B2", 'Shared region' = "#E69F00", 'knob 180' = "#D55E00", 'kindr' = "#CC79A7", 'kin10-like' =  "#009E73")

#This determines the grid colors
grid_color=rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrBr")))(500))

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
    TEMP_EXP <- DF_EXTRA_3$Name
    #These lines determine the number of iterations needed to cluster each experimental sample, while only introduces 25% experimental samples with all others being control samples. q is the rounded down whole number, r is the remainder. s_EXP is a vector repeating 25% of the number of controls q times and then adding r to 25% the number of controls for the last number. 
    q <- floor((length(TEMP_EXP)/(ncol(heat_data_control)*.25)))
    r <- length(TEMP_EXP)-(q*(ncol(heat_data_control)*.25))
    s_EXP <- c(rep((ncol(heat_data_control)*.25),q), (ncol(heat_data_control)*.25)+r)
    
    #This creates a data frame to add the Kmeans groupings too
    KMeans_Groups_All <- data.frame(Name = c(NA), Group = c(NA), Ab10_Status = c(NA), Correct = c(NA))
    KMeans_Groups_All <- KMeans_Groups_All[-c(1),]
    
    #This randomly selects the appropriate number of experimental samples as determined above, adds them to the controls, and writes out the file. It then removes the selected experimental samples from the pool and repeat's the process. In this way each experimental sample is only selected once. 
    for(i in 1:q+1) {
      #This pulls the appropriate number of experimental samples for this iteration
      n_EXP <- s_EXP[i]
      #This randomly selects n_EXP experimental samples from the vector or all experimental samples without replacement by index value. 
      samp_EXP <- sample(1:length(TEMP_EXP), n_EXP, replace = FALSE)
      #This pulls the randomly selected index values from the TEMP_EXP vector, returning actual sample names and subsets the large heat_data matrix to those columns by name. 
      heat_data_exp <- heat_data[,c(TEMP_EXP[samp_EXP])]
      #This adds the experimental samples to the control samples
      heat_data_samp <- cbind(heat_data_control, heat_data_exp)
      #This writes out the file
      write.table(heat_data_samp, file = paste("Ab10_Model/Experimental_Ab10N10Model/SubSamples/Ab10N10_Exp_heat_data_samp_MixedControls", IT, i, "txt", sep ="."), row.names = FALSE, quote = FALSE, sep = "\t")
      #This removes the experimental samples from the list of available experimental samples so that they are not chosen again. 
      TEMP_EXP <- TEMP_EXP[-samp_EXP]
    }
    
    k=1
    while(k <= q+1) {
      print(k)
      heat_data_samp <- read.table(paste("Ab10_Model/Experimental_Ab10N10Model/SubSamples/Ab10N10_Exp_heat_data_samp_MixedControls", IT, k, "txt", sep ="."), header = TRUE,  sep = "\t", check.names=FALSE)
      df_col_samp <- df_col[which(df_col$Name %in% colnames(heat_data_samp)),]
      heat_data_samp <- heat_data_samp[,df_col_samp$Name]
      #This groups by Kmeans
      group = kmeans(t(heat_data_samp[-c(13:17, 20:27),]), centers = 2)$cluster
      
      #This interprets and writes out the kmeans 
      KMeans_Groups <-  data.frame(Name=colnames(heat_data_samp), Group=group)
      
      KMeans_Groups_2 <- merge(KMeans_Groups, GROUPS, by = "Name")
      KMeans_Groups_3 <- KMeans_Groups_2[,c("Name", "Group.x", "Ab10_Status")]
      colnames(KMeans_Groups_3) <- c("Name", "Group", "Ab10_Status")
      
      KMeans_Groups_3$Ab10_Status <-  gsub("Ab10-III", "Ab10", KMeans_Groups_3$Ab10_Status)
      KMeans_Groups_3$Ab10_Status <-  gsub("Ab10-II", "Ab10", KMeans_Groups_3$Ab10_Status)
      KMeans_Groups_3$Ab10_Status <-  gsub("Ab10-I", "Ab10", KMeans_Groups_3$Ab10_Status)
      
      #This calculates the percent of each Kmeans group that is called Ab10 or N10
      Percent <- KMeans_Groups_3 %>%
        group_by(Group) %>% 
        count(Ab10_Status) %>%
        group_by(Group) %>%
        mutate(percent = n/sum(n))
      
      #This determines the predominant call for each Kmeans groups
      #This does group 1
      if(length(grep(1, Percent$Group)) > 1) {
        sub <- subset(Percent, Percent$Group == 1)
        sub2 <- subset(sub, Ab10_Status != "Unknown")
        if(nrow(sub2) == 1) {
          KMeans_Groups_3$Group <- gsub(1, sub2[1,2], KMeans_Groups_3$Group)
        } else {
          KMeans_Groups_3$Group <- gsub(1, paste("Mixed", nrow(sub2), sub2[1,2], round(sub2[1,4], digits=2), sub2[2,2], round(sub2[2,4], digits=2), sep="_"), KMeans_Groups_3$Group)
        }
      } else{
        sub <- subset(Percent, Percent$Group == 1)
        KMeans_Groups_3$Group <- gsub(1, sub[1,2], KMeans_Groups_3$Group)
      }
      
      #This does group 2
      if(length(grep(2, Percent$Group)) > 1) {
        sub <- subset(Percent, Percent$Group == 2)
        sub2 <- subset(sub, Ab10_Status != "Unknown")
        if(nrow(sub2) == 1) {
          KMeans_Groups_3$Group <- gsub(2, sub2[1,2], KMeans_Groups_3$Group)
        } else {
          KMeans_Groups_3$Group <- gsub(2, paste("Mixed", nrow(sub2), sub2[1,2], round(sub2[1,4], digits=2), sub2[2,2], round(sub2[2,4], digits=2), sep="_"), KMeans_Groups_3$Group)
        }
      } else{
        sub <- subset(Percent, Percent$Group == 2)
        KMeans_Groups_3$Group <- gsub(2, sub[1,2], KMeans_Groups_3$Group)
      }
      
      #This determines control correctness 
      KMeans_Groups_4 <- subset(KMeans_Groups_3, Ab10_Status != "Unknown")
      KMeans_Groups_4$Correct <- ifelse(KMeans_Groups_4$Group == KMeans_Groups_4$Ab10_Status, 1, 0)
      
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
      
      pdf(file = paste("Ab10_Model/Experimental_Ab10N10Model/Images/Ab10N10_HeatMap_MixedControls", IT, j, k,"pdf", sep="."))
      a <- Heatmap(as.matrix(heat_data_samp),
                   #plotting
                   col = grid_color,
                   column_title = "Tag Index Across Ab10 Haplotype",
                   column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                   row_title = "Ab10 Haplotype",
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

fwrite(Final, file = paste("Ab10_Model/Experimental_Ab10N10Model/KmeansGroups_Files/Kmeans_Exp_Groups_MixedControls", IT, "table", sep="."))

fwrite(QC_Controls, file = paste("Ab10_Model/Experimental_Ab10N10Model/QualityControl/QualityControl_KMeansGroupingControls_MixedControls", IT, "table", sep="."))

