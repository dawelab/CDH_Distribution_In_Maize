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

#Loads in the groups file
GROUPS <- vroom::vroom("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Data/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This loads in the filtered, Romero Navarro merged data set. 
B_COMB <- vroom::vroom("BChrom_TagSumPlusTagDensAllSamples.csv")

#This loads in the BINS file
BINS <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Bins_NoOverlap_BChrom.table.xlsx")

#This loads in the groups file with modified names
DF <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table")

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

#This subsets the head_data file to only the lines selected above
heat_data_control <- heat_data[,DF_EXTRA_2$Name]

#This generates the dataframe to annotate the colmns
df_col <- DF_EXTRA_2[,c("Name", "B_Chrom_Status", "Data_Source")]
colnames(df_col) <- c("Name", "B Chrom", "Data Source")
df_col$`B Chrom` <- as.factor(df_col$`B Chrom`)
df_col$`Data Source` <- as.factor(df_col$`Data Source`)

#This section reorders the column annotation by the chromosome 10 type
df_col <- df_col[order(df_col$`B Chrom`),]

#This line reorders the heat_data to match the column annotation
heat_data_control <- heat_data_control[,df_col$Name]

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

df_row <- as.data.frame(df_row$feature)
rownames(df_row) = rownames(C_bed)
colnames(df_row) <- "Feature"
df_row$Feature <- as.factor(df_row$Feature)

#This defines low copy b chromosme samples for later selection
Low <- c("NSL-2833_B-Chrom.3.DC2", "B542C_L289_B-Chrom.2.DC2", "B542C_L289_B-Chrom.3.DC2", "B542C_L289_B-Chrom.4.DC2")


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
Final <- data.frame(Name = colnames(heat_data_control))

#This defines rows that consistently have coverage in all samples
Drop_Rows <- c(4, 107, 2, 5, 87)

for(j in 1:100) {
  #This is the iteration
  print(paste("This is iteration number", j, sep=" "))
  #This ensures that roughly equal amounts of N10 and Ab10 are selected each time 
  TEMP <- merge(data.frame(Name = colnames(heat_data_control)), GROUPS)
  TEMP_No <- subset(TEMP, B_Chrom_Status == "No")
  s_No <- c(158)
  TEMP_Yes <- subset(TEMP, B_Chrom_Status == "Yes")
  #This drops low copy B Chromosomes
  TEMP_Yes <- TEMP_Yes[c(which(TEMP_Yes$Name %in% Low)),]
  s_Yes <- c(4)
  
  #This creates a data frame to add the Kmeans groupings too
  KMeans_Groups_All <- data.frame(Name = c(NA), Group = c(NA), B_Chrom_Status = c(NA), Correct = c(NA))
  KMeans_Groups_All <- KMeans_Groups_All[-c(1),]
  
  #This separates the controls into 3 groups that are roughly equal in Ab10 and N10 proportion
  for(i in 1:1) {
    n_No <- s_No[i]
    n_Yes <- s_Yes[i]
    samp_No <- sample(1:length(TEMP_No$Name), n_No, replace = FALSE)
    samp_Yes <- sample(1:length(TEMP_Yes$Name), n_Yes, replace = FALSE)
    heat_data_samp <- heat_data_control[,c(TEMP_No$Name[samp_No], TEMP_Yes$Name[samp_Yes])]
    write.table(heat_data_samp, file = paste("heat_data_samp", i, "txt", sep ="."), row.names = FALSE, quote = FALSE, sep = "\t")
    TEMP_No <- TEMP_No[-samp_No,]
  }
  
  #This groups the sub sample by kmeans, determines the Ab10/N10 identity of each group, determines the correctness of the Kmeans call, writes out the data, and plots it. 
  for(k in 1:1) {
    #This loads in the samples and subsamples the column annotations
    heat_data_samp <- read.table(paste("heat_data_samp", k, "txt", sep ="."), header = TRUE,  sep = "\t", check.names=FALSE)
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
    KMeans_Groups_3$Correct <- ifelse(KMeans_Groups_3$Group == KMeans_Groups_3$B_Chrom_Status, 1, 2)
    #This stores the calls for this run 
    KMeans_Groups_All <<- rbind(KMeans_Groups_All, KMeans_Groups_3)
    
    #This makes the heatmap annotations and plots it
    ha_col = HeatmapAnnotation(`B Chrom` = df_col_samp$`B Chrom`,
                               `Data Source` = df_col_samp$`Data Source`,
                               col = list(`B Chrom` = columnBStat_colors, `Data Source` = columndatasource_colors))
    
    ha_row = rowAnnotation(Feature = df_row$Feature,
                           col = list(Feature = row_colors),
                           na_col = "white"
    )
    
    pdf(file = paste("BChrom_Controls_HeatMap", j, k,"pdf", sep="."))
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
                 show_row_names=FALSE,
                 show_column_names = TRUE,
                 column_names_gp = grid::gpar(fontsize = 1),
    )
    print(a)
    dev.off()
  }
  #This stores the calls and writes them out
  KMeans_Groups_All_v2 <- KMeans_Groups_All
  colnames(KMeans_Groups_All_v2) <- c("Name", "Group", "B_Chrom_Status", paste("Correct", j, sep="_"))
  Final <<- merge(Final, KMeans_Groups_All_v2[,c(1,4)], by="Name")
}

write.csv(Final, file = "BChromLow_Kmeans_Groups_Final.csv", row.names = FALSE, quote = FALSE)

Correct_Sum <- Final

Correct_Sum$CorrectCount <- rowSums(Correct_Sum[,-c(1)])
Correct_Sum$Total <- ncol(Correct_Sum)-2
Correct_Sum$PropCorrect <- Correct_Sum$CorrectCount/Correct_Sum$Total
Correct_Sum_v2 <- Correct_Sum[,c(1,ncol(Correct_Sum))]
Correct_Sum_v2 <- merge(Correct_Sum_v2, GROUPS)

png("PropCorrect_BChrom_Low.png")
b <- ggplot(data=Correct_Sum_v2, aes(x=B_Chrom_Status, y=PropCorrect, color=Maize_Type)) +
  geom_jitter(width = NULL, height = 0,) +
  ylim(0,1)
print(b)
dev.off()
