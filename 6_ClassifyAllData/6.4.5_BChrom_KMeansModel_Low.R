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
setwd("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChrom_Model/Experimental_BChromModel_Low")

#This loads the original groups file
GROUPS <- vroom::vroom("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChrom_Model/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This loads in the groups file with modified names
DF <- vroom::vroom("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChrom_Model/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table")

#This loads in the BINS file
BINS <- read_excel("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChrom_Model/Bins_NoOverlap_BChrom.table.xlsx")

#This loads in the filtered, Romero Navarro merged data set. 
B_COMB <- vroom::vroom("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChrom_Model/BChrom_TagSumPlusTagDensAllSamples.csv")

BHigh <- vroom::vroom("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChrom_Model/Experimental_BChromModel_High/KmeansGroups_Files/Kmeans_Exp_Groups_High_All.csv")

BHigh_edit <- as.data.frame(apply(BHigh, 2, function(x) {gsub("Yes", 1, x)}))
BHigh_edit <- as.data.frame(apply(BHigh_edit, 2, function(x) {gsub("No", 0, x)}))
BHigh_edit2 <- as.data.frame(apply(BHigh_edit[,-c(1)], 2, function(x) {as.numeric(x)}))
BHigh_edit2$Name <- BHigh_edit$Name
BHigh_edit2$All <- rowSums(BHigh_edit2[,-c(ncol(BHigh_edit2))])
BHigh_edit2$Prop <- BHigh_edit2$All/125

test <- BHigh_edit2[,c("Name", "Prop")]

ggplot(data=BHigh_edit2, aes(x=Prop)) +
  geom_histogram()
ggsave("B_Chromosome_High_Prop_Calls.png")

ggplot(data=BHigh_edit2, aes(x=Prop)) +
  geom_histogram() +
  ylim(0,400)
ggsave("B_Chromosome_High_Prop_Calls_Zoom.png")

BHigh_final <- subset(BHigh_edit2, Prop >= 0.80)
BHigh_final <- BHigh_final[,c("Name", "Prop")]

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

#This defines low copy b chromosmes and maintains them
Low <- c("NSL-2833_B-Chrom.3.DC2", "B542C_L289_B-Chrom.2.DC2", "B542C_L289_B-Chrom.3.DC2", "B542C_L289_B-Chrom.4.DC2")
DF_EXTRA_2_No <- subset(DF_EXTRA_2, B_Chrom_Status == "No")
DF_EXTRA_2_Yes <- subset(DF_EXTRA_2, B_Chrom_Status == "Yes")
DF_EXTRA_2_Yes <- DF_EXTRA_2_Yes[c(which(DF_EXTRA_2_Yes$Name %in% Low)),]

DF_EXTRA_2 <- rbind(DF_EXTRA_2_No, DF_EXTRA_2_Yes) 

#This isolates only the experimental samples
DF_EXTRA_3 <- subset(DF_EXTRA, Data_Source != "Dawe_Lab_1" & Data_Source != "Dawe_Lab_2")

#This drops samples already called as B Chromosome in the high model
DF_EXTRA_3 <- DF_EXTRA_3[-c(which(DF_EXTRA_3$Name %in% BHigh_final$Name)),]

#This subsets the head_data file to only the lines selected above
heat_data_control <- heat_data[,DF_EXTRA_2$Name]

#This generates the dataframe to annotate the colmns
df_col <- DF_EXTRA[,c("Name", "B_Chrom_Status", "Data_Source")]
colnames(df_col) <- c("Name", "B Chrom", "Data Source")
df_col$`B Chrom` <- as.factor(df_col$`B Chrom`)
df_col$`Data Source` <- as.factor(df_col$`Data Source`)

#This generates the dataframe to annotate the rows
pos <- C_bed$start
df_row <- as.data.frame(pos)
df_row$feature <- " "

#This assigns Ab10 regions
for(i in 1:nrow(df_row)) {
  VALUE <- df_row[i,1]
  df_row[i,2] <- ifelse(VALUE >= 1 & VALUE <= 237770, 'B Short Arm', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 237771 & VALUE <= 1113960, 'B Centromere', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 1142820 & VALUE <= 8162450, 'Proximal Heterochromatin', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 8162451 & VALUE <= 20739110, 'Proximal euchromatin 1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 20739111 & VALUE <= 43049286, 'Proximal euchromatin 2', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 43049287 & VALUE <= 44287294, 'distal heterochromatin 1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 44287295 & VALUE <= 63186011, 'distal heterochromatin 2', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 63186012 & VALUE <= 71189000, 'distal heterochromatin 3', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 71189000 & VALUE <= 101843150, 'distal heterochromatin 4', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 101843151 & VALUE <= 106643106, 'distal euchromatin', df_row[i,2])
}

df_row[1,2] <- 'B Short Arm'

df_row <- as.data.frame(df_row$feature)
rownames(df_row) = rownames(C_bed)
colnames(df_row) <- "Feature"
df_row$Feature <- factor(df_row$Feature, levels= c('B Short Arm', 'B Centromere', 'Proximal Heterochromatin', 'Proximal euchromatin 1', 'Proximal euchromatin 2', 'distal heterochromatin 1', 'distal heterochromatin 2', 'distal heterochromatin 3', 'distal heterochromatin 4', 'distal euchromatin'))

#This determines the column annotation colors 
columnBStat_colors <- c("Yes" = "black", "No" ="grey80", "Unknown" = "white")

columndatasource_colors <- c("Dawe_Lab_1" = "gold", "Dawe_Lab_2" = "limegreen", "Romero-Navarro_etal_2017" = "#008080", "Swarts_etal_2017" = "#310062", "Romay_etal_2013" = "mediumorchid", "NA" = "white")

row_colors <- c('B Short Arm' = "lightpink", 'B Centromere' = "darkorange3", 'Proximal Heterochromatin' = "lightcyan2", 'Proximal euchromatin 1' = "palevioletred", 'Proximal euchromatin 2' = "hotpink", 'distal heterochromatin 1' =  "lightskyblue", "distal heterochromatin 2" = "lightskyblue3", "distal heterochromatin 3" = "deepskyblue3", "distal heterochromatin 4" = "dodgerblue4", "distal euchromatin" = "deeppink4" )

#This determines the grid colors
grid_color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500)

#heat_data_control is the matrix containing the heat data
#KMeans_Groups_All is the temporary data frame to add each run of all controls to
#j is the iteration
#Final is the data frame to add the final results of each line to 

#This creates a data frame to add the Kmeans groupings too
KMeans_Groups_All <- data.frame(Name = c(NA), Group = c(NA), Chrom_Status = c(NA), Correct = c(NA))
KMeans_Groups_All <- KMeans_Groups_All[-c(1),]

#This creates the final data frame
Final <- data.frame(Name = colnames(heat_data))

#This writes out a data drame to ensure all the controls are being classed correctly 
QC_Controls <- data.frame(Iteration=c(NA), Round=c(NA), Subround=c(NA), Controls_Correct=c(NA))
QC_Controls <- QC_Controls[-c(1),]

#This defines rows that consistently have coverage in all samples
Drop_Rows <- c(4, 107, 2, 5, 87)

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
  q <- floor((length(TEMP_EXP))/(round((ncol(heat_data_control)*.1))))
  
  r <- length(TEMP_EXP)-(q*(round((ncol(heat_data_control)*.1))))
  
  s_EXP <- c(rep(round((ncol(heat_data_control)*.1)),q-1), (round((ncol(heat_data_control)*.1))+r))
  
  
  #This creates a data frame to add the Kmeans groupings too
  KMeans_Groups_All <- data.frame(Name = c(NA), Group = c(NA), Ab10_Status = c(NA), Correct = c(NA))
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
    #write.table(heat_data_samp, file = paste("BChrom_Exp_heat_data_samp_Low", IT, i, "txt", sep ="."), row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(heat_data_samp, file = paste("SubSamples/BChrom_Exp_heat_data_samp_Low", IT, i, "txt", sep ="."), row.names = FALSE, quote = FALSE, sep = "\t")
    #This removes the experimental samples from the list of available experimental samples so that they are not chosen again. 
    TEMP_EXP <- TEMP_EXP[-samp_EXP]
  }
  
  k=1
  while(k <= q) {
    heat_data_samp <- read.table(paste("SubSamples/BChrom_Exp_heat_data_samp_Low", IT, k, "txt", sep ="."), header = TRUE,  sep = "\t", check.names=FALSE)
    #heat_data_samp <- read.table(paste("BChrom_Exp_heat_data_samp_Low", IT, k, "txt", sep ="."), header = TRUE,  sep = "\t", check.names=FALSE)
    df_col_samp <- df_col[which(df_col$Name %in% colnames(heat_data_samp)),]
    heat_data_samp <- heat_data_samp[,df_col_samp$Name]
    
    #This groups by Kmeans
    group = kmeans(t(heat_data_samp[-c(Drop_Rows),]), centers = 2)$cluster
    KMeans_Groups <-  data.frame(Name=colnames(heat_data_samp), Group=group)
    
    KMeans_Groups_2 <- merge(KMeans_Groups, GROUPS, by = "Name")
    
    KMeans_Groups_3 <- KMeans_Groups_2[,c("Name", "Group.x", "B_Chrom_Status")]
    colnames(KMeans_Groups_3) <- c("Name", "Group", "B_Chrom_Status")
    
    #This calculates the percent of each Kmeans group that is called Ab10 or N10
    Percent <- KMeans_Groups_3 %>%
      group_by(Group) %>% 
      count(B_Chrom_Status) %>%
      group_by(Group) %>%
      mutate(percent = n/sum(n))
    
    #This determines the predominant call for each Kmeans groups
    #This does group 1
    if(length(grep(1, Percent$Group)) > 1) {
      sub <- subset(Percent, Percent$Group == 1)
      sub2 <- subset(sub, B_Chrom_Status != "Unknown")
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
      sub2 <- subset(sub, B_Chrom_Status != "Unknown")
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
    KMeans_Groups_4 <- subset(KMeans_Groups_3, B_Chrom_Status != "Unknown")
    KMeans_Groups_4$Correct <- ifelse(KMeans_Groups_4$Group == KMeans_Groups_4$B_Chrom_Status, 1, 0)
    
    #This loop checks that all controls were called correctly. If they are it proceeds with the loop, if not it prints a message and redoes that same number again
    if(sum(KMeans_Groups_4$Correct) == length(KMeans_Groups_4$Correct)) {
      
      #This reports the correct controls
      line <- data.frame(Iteration=c(IT), Round=c(j), Subround=c(k), Controls_Correct=sum(KMeans_Groups_4$Correct)/nrow(KMeans_Groups_4))
      QC_Controls <<- rbind(QC_Controls, line)
      
      #This stores the calls for this run 
      KMeans_Groups_5 <- subset(KMeans_Groups_3, B_Chrom_Status == "Unknown")
      KMeans_Groups_All <<- rbind(KMeans_Groups_All, KMeans_Groups_5)
      
      #This makes the heatmap annotations
      ha_col = HeatmapAnnotation(`B Chrom` = df_col_samp$`B Chrom`,
                                 `Data Source` = df_col_samp$`Data Source`,
                                 col = list(`B Chrom` = columnBStat_colors, `Data Source` = columndatasource_colors))
      
      ha_row = rowAnnotation(Feature = df_row$Feature,
                             col = list(Feature = row_colors),
                             na_col = "white")
      
      #This rescales the data only for plotting, it does not affect grouping. Without this the differences are hard to see 
      heat_data_samp <- MinMax(heat_data_samp)
      
      #pdf(file = paste("Experimental_BChrom_HeatMap_Low", IT, j, k,"pdf", sep="."))
      pdf(file = paste("Images/BChrom_Controls_HeatMap_Low", IT, j, k,"pdf", sep="."))
      a <- Heatmap(as.matrix(heat_data_samp),
                   #plotting
                   col = grid_color,
                   column_title = "Tag Index Across B Chromosome",
                   column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                   row_title = "B Chromosome",
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
  colnames(KMeans_Groups_All_v2) <- c("Name", paste("Group", IT, j, sep="_"), "B_Chrom_Status")
  #This adds this rounds data to the output dataframe
  KMeans_Groups_All_v2_sub <- KMeans_Groups_All_v2[,c(1,2)]
  Final <<- merge(Final,KMeans_Groups_All_v2_sub, by="Name")
}

fwrite(Final, file = paste("KmeansGroups_Files/Kmeans_Exp_Groups_Low", IT, "table", sep="."))

fwrite(QC_Controls, file = paste("QualityControl/QualityControl_KMeansGroupingControls_Low", IT, "table", sep="."))
