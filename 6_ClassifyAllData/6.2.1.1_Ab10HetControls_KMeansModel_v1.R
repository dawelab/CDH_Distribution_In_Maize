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

#This sets the working directory 
setwd("")

#Loads in the groups file
#This is from 1.7
GROUPS <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This loads in the filtered, Romero Navarro merged data set.
#This is from 3.4
MERGE_Ab10Hap_RPM_FILT_2_FIX <- vroom::vroom("Ab10_Model/BWAaln_All_v_B73-Ab10_BChrom.Ab10.RPM.RNMean.table")

#This loads in the BINS file
#This is from 4.1
BINS <- read_excel("Bins_NoOverlap.table.xlsx")

#This loads in the groups file with modified names
#This is from 6.1.1
DF <- vroom::vroom("Ab10_Model/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table")

#This function goes over each row and divides each value by the max in that row 
MinMax = function(xx) { sweep(xx, 1, apply(xx, 1, max), '/') }

######## Tag Sums
#This sums the number of tag numbers across bins, but excludes the bin column from the summing

A_SUM = as.data.frame(apply(MERGE_Ab10Hap_RPM_FILT_2_FIX[,-c(1:7)], 2, function(xx) { by(xx, MERGE_Ab10Hap_RPM_FILT_2_FIX$bin, sum)  }))

#This writes out the sumed file
fwrite(A_SUM, file ="Ab10_Model/Controls_Ab10N10Model/Ab10_TagSums_AllSamples.csv")

######## Tag Density
#This sums the number of tag numbers across bins, but excludes the bin column from the summing
A_DENS = as.data.frame(apply(MERGE_Ab10Hap_RPM_FILT_2_FIX[,-c(1:7)], 2, function(xx) { by(xx, MERGE_Ab10Hap_RPM_FILT_2_FIX$bin, function(yy) { sum(yy > 0) }) }))

fwrite(A_DENS, file ="Ab10_Model/Controls_Ab10N10Model/Ab10_TagDens_AllSamples.csv")

#This takes the square root of each value
B_DENS <- sqrt(A_DENS)

#This combines the density and sum data
B_COMB <- A_SUM+B_DENS

fwrite(B_COMB, file ="Ab10_Model/Controls_Ab10N10Model/Ab10_TagSumPlusTagDensAllSamples.csv")

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
DF_EXTRA_2 <- subset(DF_EXTRA, (Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2") & Ab10_Status != "K10L2" & Ab10_Status != "Ab10_Unknown" & Name != "B73_N10.1.DC2" & Name != "B73_N10.2.DC2" & Name != "B73_N10.3.DC2" & Name != "W23_AB10-I.11.DC1"& Name != "W23_AB10-I.13.DC1" & Name != "W23_AB10-II.36.DC1" & Name != "W23_N10.14.DC1")

#This subsets the head_data file to only the lines selected above
heat_data_control <- heat_data[,DF_EXTRA_2$Name]

#This generates the dataframe to annotate the colmns
df_col <- DF_EXTRA_2[,c("Name", "Ab10_Status", "Data_Source")]
colnames(df_col) <- c("Name", "Chr10 Type", "Data Source")
df_col$`Chr10 Type` <- as.factor(df_col$`Chr10 Type`)
df_col$`Data Source` <- as.factor(df_col$`Data Source`)

#This section reorders the column annotation by the chromosome 10 type
df_col <- df_col[order(df_col$`Chr10 Type`),]

#This line reorders the heat_data to match the column annotation
heat_data_control <- heat_data_control[,df_col$Name]

#This generates the dataframe to annotate the rows
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

#This defines the homozygous lines. These need to be clustered separately. 
HOM <- c("FFMM_Ab10-I_hom.1.DC2", "FFMM_Ab10-I_hom.2.DC2", "FFMM_Ab10-I_hom.3.DC2", "W23_AB10-I.1.DC1", "W23_AB10-I.36.DC1", "W23_AB10-I.17.DC1", "W23_AB10-I.38.DC1", "W23_AB10-I.5.DC1", "W23_AB10-I.8.DC1", "W23_AB10-I.3.DC1", "W23_AB10-I.4.DC1", "W23_AB10-II.34.DC1", "W23_AB10-II.29.DC1", "W23_AB10-II.14.DC1", "W23_AB10-II.6.DC1", "W23_AB10-II.28.DC1", "W23_AB10-II.33.DC1", "W23_AB10-II.1.DC1", "W23_AB10-II.11.DC1", "W23_AB10-II.17.DC1", "W23_AB10-II.25.DC1", "W23_AB10-II.21.DC1", "W23_AB10-II.24.DC1", "W23_AB10-II.19.DC1", "W23_AB10-II.10.DC1", "W23_AB10-II.5.DC1", "W23_AB10-II.22.DC1", "W23_AB10-II.30.DC1")

#This determines the column annotation colors 
columntype_colors <- c("Ab10-I" = "black","Ab10-II" = "grey20","Ab10-III" =  "grey40", "N10" ="grey80", "Unknown" = "white")

columndatasource_colors <- c("Dawe_Lab_1" = "gold", "Dawe_Lab_2" = "limegreen", "Romero-Navarro_etal_2017" = "#008080", "Swarts_etal_2017" = "#310062", "Romay_etal_2013" = "mediumorchid", "NA" = "white")

row_colors <- c('TR1' = "#56B4E9", 'trkin' = "#0072B2", 'Shared region' = "#E69F00", 'knob 180' = "#D55E00", 'kindr' = "#CC79A7", 'kin10-like' =  "#009E73")

#This determines the grid colors
grid_color=rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrBr")))(500))



#heat_data_control is the matrix containing the heat data
#KMeans_Groups_All is the temporary data frame to add each run of all controls to
#j is the iteration
#Final is the data frame to add the final results of each line to 

#This creates a data frame to add the Kmeans groupings too
KMeans_Groups_All <- data.frame(Name = c(NA), Group = c(NA), Ab10_Status = c(NA), Correct = c(NA))
KMeans_Groups_All <- KMeans_Groups_All[-c(1),]

#This creates the final data frame
Final <- data.frame(Name = colnames(heat_data_control))


for(j in 1:100) {
  #This is the iteration
  print(paste("This is iteration number", j, sep=" "))
  #This ensures that roughly equal amounts of N10 and Ab10 are selected each time 
  TEMP <- merge(data.frame(Name = colnames(heat_data_control)), GROUPS)
  TEMP_N10 <- subset(TEMP, Ab10_Status == "N10")
  s_N10 <- c(rep(15, 3))
  TEMP_Ab10 <- subset(TEMP, Ab10_Status != "N10")
  #This drops any homozygous Ab10 lines
  TEMP_Ab10 <- TEMP_Ab10[-c(which(TEMP_Ab10$Name %in% HOM)),]
  s_Ab10 <- c(31,31, 33)
  
  #This creates a data frame to add the Kmeans groupings too
  KMeans_Groups_All <- data.frame(Name = c(NA), Group = c(NA), Ab10_Status = c(NA), Correct = c(NA))
  KMeans_Groups_All <- KMeans_Groups_All[-c(1),]
  
  #This separates the controls into 3 groups that are roughly equal in Ab10 and N10 proportion
  for(i in 1:3) {
    n_N10 <- s_N10[i]
    n_Ab10 <- s_Ab10[i]
    samp_N10 <- sample(1:length(TEMP_N10$Name), n_N10, replace = FALSE)
    samp_Ab10 <- sample(1:length(TEMP_Ab10$Name), n_Ab10, replace = FALSE)
    heat_data_samp <- heat_data_control[,c(TEMP_N10$Name[samp_N10], TEMP_Ab10$Name[samp_Ab10])]
    write.table(heat_data_samp, file = paste("Ab10_Model/Controls_Ab10N10Model/SubSamples/heat_data_samp", i, "txt", sep ="."), row.names = FALSE, quote = FALSE, sep = "\t")
    TEMP_N10 <- TEMP_N10[-samp_N10,]
    TEMP_Ab10 <- TEMP_Ab10[-samp_Ab10,]
  }
  
  #This groups the sub sample by kmeans, determines the Ab10/N10 identity of each group, determines the correctness of the Kmeans call, writes out the data, and plots it. 
  for(k in 1:3) {
    #This loads in the samples and subsamples the column annotations
    heat_data_samp <- read.table(paste("Ab10_Model/Controls_Ab10N10Model/SubSamples/heat_data_samp", i, "txt", sep ="."), header = TRUE,  sep = "\t", check.names=FALSE)
    df_col_samp <- df_col[which(df_col$Name %in% colnames(heat_data_samp)),]
    heat_data_samp <- heat_data_samp[,df_col_samp$Name]
    
    #This groups by Kmeans
    group = kmeans(t(heat_data_samp[-c(13:17, 20:27),]), centers = 2)$cluster
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
      max <- max(sub$percent)
      sub2 <- sub[grep(max, sub$percent),]
      if(sub2[1,4] >= 0.8) {
        KMeans_Groups_3$Group <- gsub(1, sub2[1,2], KMeans_Groups_3$Group)
      }
    } else {
      sub <- subset(Percent, Percent$Group == 1)
      KMeans_Groups_3$Group <- gsub(1, sub[1,2], KMeans_Groups_3$Group)
    }
    
    #This does group 2
    if(length(grep(2, Percent$Group)) > 1) {
      sub <- subset(Percent, Percent$Group == 2)
      max <- max(sub$percent)
      sub2 <- sub[grep(max, sub$percent),]
      if(sub2[1,4] >= 0.8) {
        KMeans_Groups_3$Group <- gsub(2, sub2[1,2], KMeans_Groups_3$Group)
      }
    } else {
      sub <- subset(Percent, Percent$Group == 2)
      KMeans_Groups_3$Group <- gsub(2, sub[1,2], KMeans_Groups_3$Group)
    }
    
    #This determines correctness 
    KMeans_Groups_3$Correct <- ifelse(KMeans_Groups_3$Group == KMeans_Groups_3$Ab10_Status, "Correct", "Incorrect")
    #This stores the calls for this run 
    KMeans_Groups_All <<- rbind(KMeans_Groups_All, KMeans_Groups_3)
    
    #This makes the heatmap annotations and plots it
    ha_col = HeatmapAnnotation(`Chr10 Type` = df_col_samp$`Chr10 Type`,
                               `Data Source` = df_col_samp$`Data Source`,
                               col = list(`Chr10 Type` = columntype_colors, `Data Source` = columndatasource_colors))
    
    ha_row = rowAnnotation(Feature = df_row$Feature,
                           col = list(Feature = row_colors),
                           na_col = "white"
    )
    pdf(file = paste("Ab10_Model/Controls_Ab10N10Model/Images/Ab10Het_Controls_HeatMap", j, k,"pdf", sep="."))
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
                 show_row_names=FALSE,
                 show_column_names = TRUE,
                 column_names_gp = grid::gpar(fontsize = 5),
    )
    print(a)
    dev.off()
  }
  #This stores the calls and writes them out
  KMeans_Groups_All_v2 <- KMeans_Groups_All
  colnames(KMeans_Groups_All_v2) <- c("Name", "Group", "Ab10_Status", paste("Correct", j, sep="_"))
  Final <<- merge(Final, KMeans_Groups_All_v2[,c(1,4)], by="Name")
}

write.csv(Final, file = "Ab10_Model/Controls_Ab10N10Model/Ab10HetControls_Kmeans_Groups_Final.csv", row.names = FALSE, quote = FALSE)
