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

########################################This loads and preps the data 
GROUPS <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

B_COMB <- vroom::vroom("K10L2_TagSumPlusTagDensAllSamples.csv")

#This loads in the groups file with modified names
DF <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table")

#This loads in the BINS file
BINS <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Bins_NoOverlap_K10L2.table.xlsx")

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
DF_EXTRA_2 <- subset(DF_EXTRA, (Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2") & Ab10_Status != "Ab10_Unknown" & Ab10_Status != "Ab10-II" & Name != "W23_AB10-I.11.DC1"& Name != "W23_AB10-I.13.DC1" & Name != "W23_AB10-II.36.DC1" & Name != "W23_AB10-II.37.DC1" & Name != "W23_N10.14.DC1")


#This subsets only the controls 
heat_data_control <- heat_data[,DF_EXTRA_2$Name]

#Select only the relevant columns
df_col <- DF_EXTRA_2[,c("Name", "Ab10_Status", "Data_Source")]
colnames(df_col) <- c("Name", "Chr10 Type", "Data Source")
df_col$`Chr10 Type` <- as.factor(df_col$`Chr10 Type`)
df_col$`Data Source` <- as.factor(df_col$`Data Source`)

#This section reorders the data by their chromosome 10 status
df_col <- df_col[order(df_col$`Chr10 Type`),]

#This line reorders the heat_data to match the df_col
heat_data_control <- heat_data_control[,df_col$Name]

#This creates a dataframe for the position data 
pos <- C_bed$start
df_row <- as.data.frame(pos)
df_row$feature <- " "

#This assigns Ab10 regions
for(i in 1:nrow(df_row)) {
  VALUE <- df_row[i,1]
  df_row[i,2] <- ifelse(VALUE >= 0 & VALUE <= 1381015, 'Shared region', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 1383651 & VALUE <= 7920490, 'TR1', df_row[i,2])
  df_row[i,2] <- ifelse(VALUE >= 10652268 & VALUE <= 19464951, 'TR1', df_row[i,2])
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
grid_color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(600)

#This arranges the controls for plotting
df_col_K10L2 <- df_col[which(df_col$Name %in% colnames(heat_data_control)),]
heat_data_K10L2 <- heat_data_control[,df_col_K10L2$Name]

#This makes the heatmap annotations and plots it
ha_col = HeatmapAnnotation(`Chr10 Type` = df_col_K10L2$`Chr10 Type`,
                           `Data Source` = df_col_K10L2$`Data Source`,
                           col = list(`Chr10 Type` = columntype_colors, `Data Source` = columndatasource_colors))

ha_row = rowAnnotation(Feature = df_row$Feature,
                       col = list(Feature = row_colors),
                       na_col = "white"
)

pdf(file = paste("K10L2_Controls_HeatMap_Ab10K10L2_AllWarD", "pdf", sep="."))
a <- Heatmap(as.matrix(heat_data_K10L2),
             #plotting
             col = grid_color,
             column_title = "Tag Index Across K10L2 Haplotype",
             column_title_gp = gpar(fontsize = 20, fontface = "bold"),
             row_title = "K10L2 Haplotype", 
             border_gp = gpar(col = "black", lty = 1),
             
             #clustering
             cluster_rows = FALSE, 
             cluster_columns = TRUE,
             clustering_method_columns = "ward.D",
             column_gap = unit(5, "mm"),
             
             #annotations
             top_annotation = ha_col,
             left_annotation = ha_row, 
             show_row_names=FALSE,
             show_column_names = TRUE,
             column_names_gp = grid::gpar(fontsize = 3),
)
print(a)
dev.off()


# I determined that two samples were K10L2 homozygous based on the above plots
HOM <- c("FFMM_Ab10-I_hom.1.DC2", "FFMM_Ab10-I_hom.2.DC2", "FFMM_Ab10-I_hom.3.DC2", "W23_AB10-I.1.DC1", "W23_AB10-I.36.DC1", "W23_AB10-I.17.DC1", "W23_AB10-I.38.DC1", "W23_AB10-I.5.DC1", "W23_AB10-I.8.DC1", "W23_AB10-I.3.DC1", "W23_AB10-I.4.DC1", "W23_AB10-II.34.DC1", "W23_AB10-II.29.DC1", "W23_AB10-II.14.DC1", "W23_AB10-II.6.DC1", "W23_AB10-II.28.DC1", "W23_AB10-II.33.DC1", "W23_AB10-II.1.DC1", "W23_AB10-II.11.DC1", "W23_AB10-II.17.DC1", "W23_AB10-II.25.DC1", "W23_AB10-II.21.DC1", "W23_AB10-II.24.DC1", "W23_AB10-II.19.DC1", "W23_AB10-II.10.DC1", "W23_AB10-II.5.DC1", "W23_AB10-II.22.DC1", "W23_AB10-II.30.DC1", "B73_K10L2.10.DC2", "B73_K10L2.4.DC2")

heat_data_K10L2_het <- heat_data_K10L2[,-c(which(colnames(heat_data_K10L2) %in% HOM))]
df_col_K10L2_het <- df_col[which(df_col$Name %in% colnames(heat_data_K10L2_het)),]
heat_data_K10L2_het <- heat_data_control[,df_col_K10L2_het$Name]

set.seed(123)
group = kmeans(t(heat_data_K10L2_het[c(4:15),]), centers = 2)$cluster

#This makes the heatmap annotations and plots it
ha_col = HeatmapAnnotation(`Chr10 Type` = df_col_K10L2_het$`Chr10 Type`,
                           `Data Source` = df_col_K10L2_het$`Data Source`,
                           col = list(`Chr10 Type` = columntype_colors, `Data Source` = columndatasource_colors))



pdf(file = paste("K10L2_Controls_HeatMap_Ab10K10L2_AllKMeans", "pdf", sep="."))
a <- Heatmap(as.matrix(heat_data_K10L2_het),
             #plotting
             col = grid_color,
             column_title = "Tag Index Across K10L2 Haplotype",
             column_title_gp = gpar(fontsize = 20, fontface = "bold"),
             row_title = "K10L2 Haplotype", 
             border_gp = gpar(col = "black", lty = 1),
             
             #clustering
             cluster_rows = FALSE, 
             cluster_columns = TRUE,
             column_split = group,
             column_gap = unit(5, "mm"),
             
             #annotations
             top_annotation = ha_col,
             left_annotation = ha_row, 
             show_row_names=FALSE,
             show_column_names = TRUE,
             column_names_gp = grid::gpar(fontsize = 3),
)
print(a)
dev.off()


#This creates a data frame to add the Kmeans groupings too
KMeans_Groups_All <- data.frame(Name = c(NA), Group = c(NA), Ab10_Status = c(NA), Correct = c(NA))
KMeans_Groups_All <- KMeans_Groups_All[-c(1),]

#This creates the final data frame
Final_2 <- data.frame(Name = colnames(heat_data_control))

#heat_data_control is the matrix containing the heat data
#KMeans_Groups_All is the temporary data frame to add each run of all controls to
#j is the iteration
#Correct is the data frame to add the Correct results of each line to 

for(j in 1:100) {
  #This is the iteration
  print(paste("This is iteration number", j, sep=" "))
  #This ensures that roughly equal amounts of N10 and Ab10 are selected each time 
  TEMP <- merge(data.frame(Name = colnames(heat_data_control)), GROUPS)
  TEMP_N10 <- subset(TEMP, Ab10_Status == "N10")
  s_N10 <- c(24, 24)
  TEMP_K10L2 <- subset(TEMP, Ab10_Status == "K10L2")
  #This drops the homozygous samples
  TEMP_K10L2 <- TEMP_K10L2[c(which(TEMP_K10L2$Name %in% HOM)),]
  s_K10L2 <- c(1,1)
  TEMP_Ab10 <- subset(TEMP, Ab10_Status != "N10" & Ab10_Status != "K10L2" & Ab10_Status != "Ab10-II")
  #This drops the homozygous samples
  TEMP_Ab10 <- TEMP_Ab10[c(which(TEMP_Ab10$Name %in% HOM)),]
  s_Ab10 <- c(5, 6)
  
  #This creates a data frame to add the Kmeans groupings too
  KMeans_Groups_All <- data.frame(Name = c(NA), Group = c(NA), Ab10_Status = c(NA), K10L2_Status = c(NA), Correct = c(NA))
  KMeans_Groups_All <- KMeans_Groups_All[-c(1),]
  
  #This separates the controls into 3 groups that are roughly equal in Ab10 and N10 proportion
  for(i in 1:2) {
    n_N10 <- s_N10[i]
    n_K10L2 <- s_K10L2[i]
    n_Ab10 <- s_Ab10[i]
    samp_N10 <- sample(1:length(TEMP_N10$Name), n_N10, replace = FALSE)
    samp_K10L2 <- sample(1:length(TEMP_K10L2$Name), n_K10L2, replace = FALSE)
    samp_Ab10 <- sample(1:length(TEMP_Ab10$Name), n_Ab10, replace = FALSE)
    heat_data_samp <- heat_data_control[,c(TEMP_N10$Name[samp_N10], TEMP_K10L2$Name[samp_K10L2], TEMP_Ab10$Name[samp_Ab10])]
    write.table(heat_data_samp, file = paste("K10L2_heat_data_samp", i, "txt", sep ="."), row.names = FALSE, quote = FALSE, sep = "\t")
    TEMP_N10 <- TEMP_N10[-samp_N10,]
    TEMP_K10L2 <- TEMP_K10L2[-samp_K10L2,]
    TEMP_Ab10 <- TEMP_Ab10[-samp_Ab10,]
  }
  
  #This groups the sub sample by kmeans, determines the Ab10/N10 identity of each group, determines the correctness of the Kmeans call, writes out the data, and plots it. 
  
  for(k in 1:2) {
    #This loads in the samples and subsamples the column annotations
    heat_data_samp <- read.table(paste("K10L2_heat_data_samp", k, "txt", sep ="."), header = TRUE,  sep = "\t", check.names=FALSE)
    df_col_samp <- df_col[which(df_col$Name %in% colnames(heat_data_samp)),]
    heat_data_samp <- heat_data_samp[,df_col_samp$Name]
    
    #This groups by Kmeans
    group = kmeans(t(heat_data_samp[c(7:10),]), centers = 2)$cluster
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
      max <- max(sub$percent)
      sub2 <- sub[grep(max, sub$percent),]
      if(sub2[1,4] >= 0.8) {
        KMeans_Groups_3$Group <- gsub("\\b1\\b", sub2[1,2], KMeans_Groups_3$Group, perl = TRUE)
      }
    } else {
      sub <- subset(Percent, Percent$Group == 1)
      KMeans_Groups_3$Group <- gsub("\\b1\\b", sub[1,2], KMeans_Groups_3$Group, perl = TRUE)
    }
    
    #This does group 2
    if(length(grep(2, Percent$Group)) > 1) {
      sub <- subset(Percent, Percent$Group == 2)
      max <- max(sub$percent)
      sub2 <- sub[grep(max, sub$percent),]
      if(sub2[1,4] >= 0.8) {
        KMeans_Groups_3$Group <- gsub("\\b2\\b", sub2[1,2], KMeans_Groups_3$Group, perl = TRUE)
      }
    } else {
      sub <- subset(Percent, Percent$Group == 2)
      KMeans_Groups_3$Group <- gsub("\\b2\\b", sub[1,2], KMeans_Groups_3$Group, perl = TRUE)
    }
    
    #This determines correctness 
    KMeans_Groups_3$Correct <- ifelse(KMeans_Groups_3$Group == KMeans_Groups_3$K10L2_Status, 1, 0)
    
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
    pdf(file = paste("K10L2_Controls_HeatMap", j, k,"pdf", sep="."))
    a <- Heatmap(as.matrix(heat_data_samp),
                 #plotting
                 col = grid_color,
                 column_title = "Tag Index Across Ab10 Haplotype",
                 column_title_gp = gpar(fontsize = 20, fontface = "bold"),
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
                 show_row_names=FALSE,
                 show_column_names = TRUE,
                 column_names_gp = grid::gpar(fontsize = 5),
    )
    print(a)
    dev.off()
  }
  
  #This stores the calls and writes them out
  KMeans_Groups_All_v2 <- KMeans_Groups_All
  colnames(KMeans_Groups_All_v2) <- c("Name", "Group", "Ab10_Status", "K10L2_Status", paste("Correct", j, sep="_"))
  Final_2 <<- merge(Final_2, KMeans_Groups_All_v2[,c(1,ncol(KMeans_Groups_All_v2))], by="Name")
}

Correct_Sum <- Final_2

Correct_Sum$CorrectCount <- rowSums(Correct_Sum[,-c(1)])
Correct_Sum$Total <- ncol(Correct_Sum)-2
Correct_Sum$PropCorrect <- Correct_Sum$CorrectCount/Correct_Sum$Total
Correct_Sum_v2 <- Correct_Sum[,c(1,ncol(Correct_Sum))]
Correct_Sum_v2 <- merge(Correct_Sum_v2, GROUPS)

png("PropCorrect_K10L2Ab10_Hom_7-10.png")
b <- ggplot(data=Correct_Sum_v2, aes(x=Ab10_Status, y=PropCorrect, color=Maize_Type)) +
  geom_jitter(width = NULL, height = 0,) +
  ylim(0,1)
print(b)
dev.off()
