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
library(randomForest)
library(caret)
library(NbClust)


########################################This loads and preps the data 
#This is from 1.7
GROUPS <- vroom::vroom("Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

#This is from 3.4
MERGE_Ab10Hap_RPM_FILT_2_FIX <- vroom::vroom("BWAaln_All_v_Ab10HIFIBChrom.Ab10.RPM.RNMean.table")

#This is from 6.1.1
DF <- vroom::vroom("Ab10Model/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_NameChanges.table")

#This loads in the BINS file
#This is from 4.1
BINS <- read_excel("Ab10_Global_Survey/Bins_NoOverlap.table.xlsx")

#This is from 6.2.3
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
DF_EXTRA_2 <- subset(DF_EXTRA, (Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2") & Ab10_Status != "K10L2" & Ab10_Status != "N10" & Ab10_Status != "Ab10_Unknown" & Name != "B73_N10.1.DC2" & Name != "B73_N10.2.DC2" & Name != "B73_N10.3.DC2" & Name != "W23_AB10-I.11.DC1"& Name != "W23_AB10-I.13.DC1" & Name != "W23_AB10-II.36.DC1" & Name != "W23_N10.14.DC1" & Name != "Ab10-I-Jal_het.1.DC2" & Name != "Ab10-I-Jal_het.2.DC2" & Name != "Ab10-I-Jal_het.3.DC2" & Name != "Ab10-I-Jal_het.4.DC2" & Name != "Ab10-I-Jal_het.5.DC2")

#This subsets only the Ab10 controls 
heat_data_control <- heat_data[,DF_EXTRA_2$Name]

#I can separate Ab10I and III from Ab10 II perfectly using the whole haplotype. I can't separate Ab10-I and III 100% there are always a few wrong

#This transposes the heat data and merges it with the Ab10 status

forRF_temp <- as.data.frame(t(heat_data_control))
forRF_temp$Name <- rownames(forRF_temp)
forRF_temp <- merge(forRF_temp, GROUPS[,c("Name", "Ab10_Status")], )
rownames(forRF_temp) <- forRF_temp$Name
forRF_temp <- forRF_temp[,-c(1)]
forRF <- forRF_temp
rownames(forRF) <- rownames(forRF_temp)

#This selects Ab10 positive experimental samples
POS <- POS[,-c(127)]
POS$All <- rowSums(POS[,-c(1)])
POS_sub <- subset(POS, All >= 119)
Ab10_POS <- POS_sub$Name

#This splits the data into a training and test set and renames the columns as the random forest cannot have coumn names that start with numbers
  j=1
  print(j)
  set.seed(j)
  sample <- sample(c(TRUE, FALSE), nrow(forRF), replace=TRUE, prob=c(0.7,0.3))
  train <- as.data.frame(forRF[sample,])
  test <- as.data.frame(forRF[!sample,])
  test$Ab10_Status <- as.factor(test$Ab10_Status)
  train$Ab10_Status <- as.factor(train$Ab10_Status)
  colnames(train) <- c("Bin_1", "Bin_2", "Bin_3", "Bin_4", "Bin_5", "Bin_6", "Bin_7", "Bin_8", "Bin_9", "Bin_10", "Bin_11", "Bin_12", "Bin_13", "Bin_14", "Bin_15", "Bin_16", "Bin_17", "Bin_18", "Bin_19", "Bin_20", "Bin_21", "Bin_22", "Bin_23", "Bin_24", "Bin_25", "Bin_26", "Bin_27", "Bin_28", "Bin_29", "Bin_30", "Bin_31", "Bin_32", "Bin_33", "Bin_34", "Bin_35", "Bin_36", "Bin_37", "Bin_38", "Bin_39", "Bin_40", "Bin_41", "Bin_42", "Bin_43", "Bin_44", "Bin_45", "Bin_46", "Bin_47", "Bin_48", "Bin_49", "Bin_50", "Bin_51", "Bin_52", "Bin_53", "Bin_54", "Ab10_Status")
 
   colnames(test) <- c("Bin_1", "Bin_2", "Bin_3", "Bin_4", "Bin_5", "Bin_6", "Bin_7", "Bin_8", "Bin_9", "Bin_10", "Bin_11", "Bin_12", "Bin_13", "Bin_14", "Bin_15", "Bin_16", "Bin_17", "Bin_18", "Bin_19", "Bin_20", "Bin_21", "Bin_22", "Bin_23", "Bin_24", "Bin_25", "Bin_26", "Bin_27", "Bin_28", "Bin_29", "Bin_30", "Bin_31", "Bin_32", "Bin_33", "Bin_34", "Bin_35", "Bin_36", "Bin_37", "Bin_38", "Bin_39", "Bin_40", "Bin_41", "Bin_42", "Bin_43", "Bin_44", "Bin_45", "Bin_46", "Bin_47", "Bin_48", "Bin_49", "Bin_50", "Bin_51", "Bin_52", "Bin_53", "Bin_54", "Ab10_Status")
  
   #This data set is unbalanced with about 30 Ab10-I and Ab10-II, but only 15 Ab10-III
  #This subset the heat data to the Ab10
  heat_data_Ab10 <- heat_data[,which(colnames(heat_data) %in% c(Ab10_POS, "Ab10-I-Jal_het.1.DC2", "Ab10-I-Jal_het.2.DC2", "Ab10-I-Jal_het.3.DC2", "Ab10-I-Jal_het.4.DC2", "Ab10-I-Jal_het.5.DC2"))]
  
  RF_Ab10 <- as.data.frame(t(heat_data_Ab10))
  
  colnames(RF_Ab10) <- c("Bin_1", "Bin_2", "Bin_3", "Bin_4", "Bin_5", "Bin_6", "Bin_7", "Bin_8", "Bin_9", "Bin_10", "Bin_11", "Bin_12", "Bin_13", "Bin_14", "Bin_15", "Bin_16", "Bin_17", "Bin_18", "Bin_19", "Bin_20", "Bin_21", "Bin_22", "Bin_23", "Bin_24", "Bin_25", "Bin_26", "Bin_27", "Bin_28", "Bin_29", "Bin_30", "Bin_31", "Bin_32", "Bin_33", "Bin_34", "Bin_35", "Bin_36", "Bin_37", "Bin_38", "Bin_39", "Bin_40", "Bin_41", "Bin_42", "Bin_43", "Bin_44", "Bin_45", "Bin_46", "Bin_47", "Bin_48", "Bin_49", "Bin_50", "Bin_51", "Bin_52", "Bin_53", "Bin_54")
  
  
  #This trains the random forest model
  rf <- randomForest(Ab10_Status~., data=train, ntree=1000000, proximity=TRUE) 
 
  #This predicts on the training data set
   p1 <- predict(rf, train)
   confusionMatrix(p1, train$Ab10_Status)
  
   #This predicts on the test data set
   p2 <- predict(rf, test)
   confusionMatrix(p2, test$Ab10_Status)
  
   #This creates several plots describing the training of the random forest
   plot(rf)
  
   hist(treesize(rf),
       main = "No. of Nodes for the Trees",
       col = "green")
  
  varImpPlot(rf,
             sort = T,
             n.var = 54,
             main = "Variable Importance")
  
  
  #This creates a dataframe for the position data 
  pos <- C_bed$start
  df_row <- as.data.frame(pos)
  df_row$feature <- " "
  
  #This creates a data frame to annotate the rows which will be used in the plotting 
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
    #df_row[i,2] <- ifelse(VALUE >= 190307761 & VALUE <= 191107857, 'kin10-like', df_row[i,2])
  }
  
  df_row <- as.data.frame(df_row$feature)
  rownames(df_row) = rownames(C_bed)
  colnames(df_row) <- "Feature"
  df_row$Feature <- as.factor(df_row$Feature)
  
  
  #This extracts information on variable importance from the random forest model 
  IMPO <- as.data.frame(importance(rf))
  IMPO <- cbind(IMPO, df_row$Feature)
  IMPO$bin <- rownames(IMPO) 
  IMPO$bin <- gsub("Bin_", "", IMPO$bin)
  
  MDSplot(rf, train$Ab10_Status)
  
  #This extracts the proportion of decision trees that called each class for each sample
  p4 <- as.data.frame(predict(rf, RF_Ab10, type="prob"))
  p4$max <- apply(p4, 1, function(x) max(x))
  p4$Name <- rownames(p4)

  #This plots the most frequent calls proportion of decision trees for all samples as a histogram
  ggplot(p4, aes(x=max)) +
    geom_histogram(bins=24) +
    labs(x="Proportion of Decision Trees With Majority Call", y= "Count") +
    scale_x_continuous(breaks=c(0.33, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
    ggtitle("Histogram of Random Forest Decision Trees With Majority Call")
ggsave("Histogram_RFPropVotesForMax.png")
  

  #This isolates samples that had more that 65% of decision trees call the same class and assigns that class
  Sub_High <- subset(p4, p4$`Ab10-I` >= 0.65 | p4$`Ab10-II` >= 0.65 | p4$`Ab10-III` >= 0.65)
  Sub_Highsmall <- Sub_High[,c("Ab10-I", "Ab10-II", "Ab10-III")]
  Sub_Highsmall$Conf_ID <- c(NA)

  for(i in 1:nrow(Sub_Highsmall)){
    index <- grep(max(Sub_Highsmall[i,c(1:3)]), Sub_Highsmall[i,c(1:3)])
    ID <- colnames(Sub_Highsmall)[index]
    Sub_Highsmall[i,4] <- ID 
  }
  
  #This isolates samples that did not have any one class called by 65% of decision trees and calls them ambiguous
  Sub_Low <- subset(p4, p4$`Ab10-I` < 0.65 & p4$`Ab10-II` < 0.65 & p4$`Ab10-III` < 0.65)
  Sub_Lowsmall <- Sub_Low[,c("Ab10-I", "Ab10-II", "Ab10-III")]
  Sub_Lowsmall$Conf_ID <- c("Ambiguous")
  
  #This merges my High and low confidence Ab10s 
  Sub_All <- rbind(Sub_Highsmall, Sub_Lowsmall)

 #This writes out the type calls 
  WRITE <-  Sub_All
  WRITE$Name <- rownames(WRITE)
  write.csv(WRITE[,c("Name", "Conf_ID")], file="RF_Ab10TypeClassifications.csv", row.names = FALSE, quote = FALSE)
                  
  #This subsets the heat data to the high confidence type calls 
  heat_data_Ab10Type <- heat_data[,rownames(Sub_All)]

  
  #Select only the relevant columns
  df_col <- DF_EXTRA[,c("Name", "Ab10_Status", "Data_Source")]
  colnames(df_col) <- c("Name", "Chr10 Type", "Data Source")
  df_col$`Chr10 Type` <- as.factor(df_col$`Chr10 Type`)
  df_col$`Data Source` <- as.factor(df_col$`Data Source`)

  #This creates the df_col dataframe 
  df_col_control <- df_col[which(df_col$Name %in% colnames(heat_data_control)),]
  
  df_col_Ab10Type <- df_col[which(df_col$Name %in% colnames(heat_data_Ab10Type)),]
  Sub_All_v2 <- Sub_All
  Sub_All_v2$Name <- rownames( Sub_All_v2)
  Sub_All_v2 <-  Sub_All_v2[, c("Name", "Conf_ID")]
  df_col_Ab10Type <- merge(df_col_Ab10Type, Sub_All_v2, by="Name")
  
  #This section adds the confidence score to the column annotation
  Sub_All_v2 <- Sub_All
  Sub_All_v2$max <- apply(Sub_All_v2[,c(1:3)], 1, function(x) max(x))
  Sub_All_v2$Name <- rownames(Sub_All_v2)
  df_col_Ab10Type <- merge(df_col_Ab10Type, Sub_All_v2, by="Name", all.x = TRUE)
  df_col_Ab10Type <- df_col_Ab10Type[,c("Name", "Chr10 Type", "Data Source", "Conf_ID.x", "max")]
  colnames(df_col_Ab10Type) <- c("Name", "Chr10 Type", "Data Source", "RF Class", "Confidence")
  
  
  #This section reorders the data to match the heat data
  df_col_control <- df_col_control[order(df_col_control$`Chr10 Type`),]
  heat_data_control <- heat_data_control[,c(df_col_control$Name)]
  
  df_col_Ab10Type <- df_col_Ab10Type[order(df_col_Ab10Type$`RF Class`),]
  heat_data_Ab10Type <- heat_data_Ab10Type[,c(df_col_Ab10Type$Name)]
  
  
  

  #This adds the weight of each variable going into the LDA model
  df_row$MeanDecreaseGini <- IMPO$MeanDecreaseGini

  #This determines the column annotation colors 
  #This determines the grid colors
  grid_color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500)

  columntypetrue_colors <- c("Ab10-I" = "#440154FF","Ab10-II" = "#2A788EFF","Ab10-III" =  "#7AD151FF", "N10" ="grey70", "Unknown" = "white")

  #columntypeLDA_colors <- c("Ab10-II" = "navy","Ab10-I" = "cornflowerblue", "Ab10-III" = "lightblue", "NA" = "seagreen")

  rowfeature_colors <- c('TR-1' = "#56B4E9", 'Trkin' = "#0072B2", 'Shared Region' = "#E69F00", 'knob180' = "#D55E00", 'Kindr' = "#CC79A7")

  rowweight_colors = colorRamp2(c(min(df_row$MeanDecreaseGini),max(df_row$MeanDecreaseGini)), c("thistle1", "darkred"))

  #This makes the heatmap annotations
  ha_col = HeatmapAnnotation(`Chr10 Type` = df_col_control$`Chr10 Type`,
                           col = list(`Chr10 Type` = columntypetrue_colors))

  ha_row = rowAnnotation(Feature = df_row$Feature,
                       `MeanDecreaseGini` = df_row$MeanDecreaseGini,
                       col = list(Feature = rowfeature_colors, `MeanDecreaseGini` = rowweight_colors),
                       na_col = "white")
  
  pdf(file = paste("Controls_HeatMap_RandomForest", "pdf", sep="."), height=10, width=10)
    
  a <- Heatmap(as.matrix(heat_data_control),
             #plotting
             col = grid_color,
             column_title = "Tag Index Across Ab10 Haplotype",
             column_title_gp = gpar(fontsize = 20, fontface = "bold"),
             row_title = "Ab10 Haplotype", 
             border_gp = gpar(col = "black", lty = 1),
             
             #clustering
             cluster_rows = FALSE, 
             cluster_columns = FALSE, 
             column_split = df_col_control$`Chr10 Type`,
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
  

  ###########This plots the Ab10 calls 
  
  #This determines the column annotation colors 
  #This determines the grid colors
  grid_color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500)
  
  columntypeRF_colors <- c("Ab10-I" = "#440154FF","Ab10-II" = "#2A788EFF","Ab10-III" =  "#7AD151FF", "Ambiguous" = "#FDE725FF")
  
  columnconf_colors = colorRamp2(c(0.65,1), c("aliceblue", "darkcyan"))
  
  rowfeature_colors <- c('TR-1' = "#56B4E9", 'Trkin' = "#0072B2", 'Shared Region' = "#E69F00", 'knob180' = "#D55E00", 'Kindr' = "#CC79A7")
  
  rowweight_colors = colorRamp2(c(min(df_row$MeanDecreaseGini),max(df_row$MeanDecreaseGini)), c("thistle1", "darkred"))
  
  #This makes the heatmap annotations
  ha_col = HeatmapAnnotation(`RF Class` = df_col_Ab10Type$`RF Class`,
                             `RF Confidence` = df_col_Ab10Type$Confidence,
                             col = list(`RF Class` =  columntypeRF_colors, `RF Confidence` = columnconf_colors), 
                             annotation_legend_param = list( `RF Confidence` = list(at=c(0.5, 0.6, 0.7, 0.8, 0.9, 1))),
                             na_col = "white"
                             )
  
  ha_row = rowAnnotation(Feature = df_row$Feature,
                         `MeanDecreaseGini` = df_row$MeanDecreaseGini,
                         col = list(Feature = rowfeature_colors, `MeanDecreaseGini` = rowweight_colors),
                         na_col = "white")
  
  pdf(file = paste("Ab10Experimental_HeatMap_RandomForest", "pdf", sep="."), height=10, width=10)
  
  a <- Heatmap(as.matrix(heat_data_Ab10Type),
               #plotting
               col = grid_color,
               column_title = paste("Tag Index Across Ab10 Haplotype\n", nrow(Sub_Highsmall), "samples classable" ),
               column_title_gp = gpar(fontsize = 20, fontface = "bold"),
               row_title = "Ab10 Haplotype", 
               border_gp = gpar(col = "black", lty = 1),
               
               #clustering
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               column_split = factor(df_col_Ab10Type$`RF Class`, level=c("Ab10-I", "Ab10-II", "Ab10-III", "Ambiguous")),
               cluster_row_slices = FALSE, 
               cluster_column_slices = FALSE,
               column_gap = unit(5, "mm"),
               show_column_dend = TRUE,
               
               
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
  
  pdf(file = paste("Ab10Experimental_HeatMap_RandomForest_colclust", "pdf", sep="."), height=10, width=10)
  
  a <- Heatmap(as.matrix(heat_data_Ab10Type),
               #plotting
               col = grid_color,
               column_title = paste("Tag Index Across Ab10 Haplotype\n", nrow(Sub_Highsmall), "samples classable" ),
               column_title_gp = gpar(fontsize = 20, fontface = "bold"),
               row_title = "Ab10 Haplotype", 
               border_gp = gpar(col = "black", lty = 1),
               
               #clustering
               cluster_rows = FALSE, 
               cluster_columns = TRUE, 
               column_split = factor(df_col_Ab10Type$`RF Class`, level=c("Ab10-I", "Ab10-II", "Ab10-III", "Ambiguous")),
               cluster_row_slices = FALSE, 
               cluster_column_slices = TRUE,
               clustering_method_column = "ward.D",
               column_gap = unit(5, "mm"),
               show_column_dend = TRUE,
               
               
               #annotations
               name="Tag Index",
               top_annotation = ha_col,
               left_annotation = ha_row, 
               show_row_names= FALSE,
               show_column_names = FALSE, 
               row_names_gp = grid::gpar(fontsize = 5),
               column_names_gp = grid::gpar(fontsize = 0.5)
  )
  print(a)
  dev.off()
  
  
  ###########This plots only the ambigious calls 
  df_col_Ab10TypeAmbig <- subset(df_col_Ab10Type, `RF Class` == "Ambiguous")
  #This pulls the closest class
  df_col_Ab10TypeAmbig <- merge(df_col_Ab10TypeAmbig, Sub_Low, by="Name")
  
  df_col_Ab10TypeAmbig$`Approximate Class` <- c(NA)
  
  i=1
  for(i in 1:nrow(df_col_Ab10TypeAmbig)){
    index <- grep(max(df_col_Ab10TypeAmbig[i,c(6:8)]), df_col_Ab10TypeAmbig[i,c(6:8)])
    ID <- colnames(df_col_Ab10TypeAmbig)[index+5]
    df_col_Ab10TypeAmbig[i,10] <- ID 
  }
  
  df_col_Ab10TypeAmbig <- df_col_Ab10TypeAmbig[,c(1:5,10)]
  
  heat_data_Ab10TypeAmbig <- heat_data_Ab10[,c(which(colnames(heat_data_Ab10) %in% df_col_Ab10TypeAmbig $Name))]
  
  #This determines the column annotation colors 
  #This determines the grid colors
  grid_color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(500)
  
  columntypetrue_colors <- c("Ab10-I" = "black","Ab10-II" = "grey20","Ab10-III" =  "grey40", "N10" ="grey80", "Unknown" = "white", "N10" = "beige")
  
  columntypeRF_colors <- c("Ab10-II" = "navy","Ab10-I" = "cornflowerblue", "Ab10-III" = "lightblue")
  
  columndatasource_colors <- c("Dawe_Lab_1" = "gold", "Dawe_Lab_2" = "limegreen", "Romero-Navarro_etal_2017" = "#008080", "Swarts_etal_2017" = "#310062", "Romay_etal_2013" = "mediumorchid", "NA" = "white")
  
  columnconf_colors = colorRamp2(c(0.3,0.65), c("aliceblue", "darkcyan"))
  
  rowfeature_colors <- c('TR1' = "#56B4E9", 'trkin' = "#0072B2", 'Shared region' = "#E69F00", 'knob 180' = "#D55E00", 'kindr' = "#CC79A7", 'kin10-like' =  "#009E73")
  
  rowweight_colors = colorRamp2(c(min(df_row$MeanDecreaseGini),max(df_row$MeanDecreaseGini)), c("thistle1", "darkred"))
  
  ha_col = HeatmapAnnotation(`Chr10 Type` = df_col_Ab10TypeAmbig$`Chr10 Type`,
                             `Data Source` = df_col_Ab10TypeAmbig$`Data Source`,
                             `Approximate RF Class` = df_col_Ab10TypeAmbig$`Approximate Class`,
                             `RF Confidence` = df_col_Ab10TypeAmbig$Confidence,
                             col = list(`Chr10 Type` = columntypetrue_colors, `Data Source` = columndatasource_colors, `Approximate RF Class` =  columntypeRF_colors, `RF Confidence` = columnconf_colors), 
                             annotation_legend_param = list( `RF Confidence` = list(at=c(0.3, 0.4, 0.5, 0.65))),
                             na_col = "white"
  )
  
  ha_row = rowAnnotation(Feature = df_row$Feature,
                         `MeanDecreaseGini` = df_row$MeanDecreaseGini,
                         col = list(Feature = rowfeature_colors, `MeanDecreaseGini` = rowweight_colors),
                         na_col = "white")
  
  pdf(file = paste("Ab10Experimental_HeatMap_RandomForest_Ambig", "pdf", sep="."), height=10, width=10)
  a <- Heatmap(as.matrix(heat_data_Ab10TypeAmbig),
               #plotting
               col = grid_color,
               column_title = "Tag Index Across Ab10 Haplotype\n Ambiguous Class",
               column_title_gp = gpar(fontsize = 20, fontface = "bold"),
               row_title = "Ab10 Haplotype", 
               border_gp = gpar(col = "black", lty = 1),
               
               #clustering
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               column_split = df_col_Ab10TypeAmbig$`RF Class`,
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
  
  pdf(file = paste("Ab10Experimental_HeatMap_RandomForest_Ambig_colclust", "pdf", sep="."), height=10, width=10)
  a <- Heatmap(as.matrix(heat_data_Ab10TypeAmbig),
               #plotting
               col = grid_color,
               column_title = "Tag Index Across Ab10 Haplotype",
               column_title_gp = gpar(fontsize = 20, fontface = "bold"),
               row_title = "Ab10 Haplotype", 
               border_gp = gpar(col = "black", lty = 1),
               
               #clustering
               cluster_rows = FALSE, 
               cluster_columns = TRUE, 
               clustering_distance_columns = "manhattan",
               column_split = df_col_Ab10TypeAmbig$`RF Class`,
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
  
  pdf(file = paste("Ab10Experimental_HeatMap_RandomForest_Ambig_colclust_rowclust", "pdf", sep="."), height=10, width=10)
  a <- Heatmap(as.matrix(heat_data_Ab10TypeAmbig),
               #plotting
               col = grid_color,
               column_title = "Tag Index Across Ab10 Haplotype",
               column_title_gp = gpar(fontsize = 20, fontface = "bold"),
               row_title = "Ab10 Haplotype", 
               border_gp = gpar(col = "black", lty = 1),
               
               #clustering
               cluster_rows = TRUE, 
               cluster_columns = TRUE, 
               clustering_distance_columns = "manhattan",
               clustering_distance_rows = "manhattan",
               column_split = df_col_Ab10TypeAmbig$`RF Class`,
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
  

  #This creates the distance matrix and then performs the PCA
  
  #This selects only landrace controls and removes Jal from the control set
  Ab10_Control_Landrace <- subset(GROUPS, Data_Source == "Dawe_Lab_2" & Maize_Type == "Landrace" & Ab10_Status != "N10" & Ab10_Status != "K10L2" & Ab10_Status != "Ab10_Unknown" & Name != "Ab10-I-Jal_het.1.DC2" & Name != "Ab10-I-Jal_het.2.DC2" & Name != "Ab10-I-Jal_het.3.DC2" & Name != "Ab10-I-Jal_het.4.DC2" & Name != "Ab10-I-Jal_het.5.DC2")
  
  heat_data_Ab10_Control_Landrace <- heat_data[,which(colnames(heat_data) %in% Ab10_Control_Landrace$Name)]
  
  heat_data_Ab10TypeAmbig_controls <- cbind(heat_data_Ab10Type, heat_data_Ab10_Control_Landrace)

  #This selects only variables from the random forest with a High mean decreasing Gini (2)
  HiGini <- c(5, 8, 11, 15, 26, 35, 36, 37, 38)
  forPCA <- t(heat_data_Ab10TypeAmbig_controls[HiGini,])
  
  #This performs the pca
  dist <- dist(forPCA)
  pca <- cmdscale(dist, k=10, eig = TRUE)
  
  #This calculates the percent of the variation they account for 
  TOT <- sum(pca$eig)
  pca$perc_var <- (pca$eig/TOT)*100 
  
  #This plots the eigen vector values. This tells you how much variability is accounted for by each PCA
  pdf("PCA_ScreePlot_Ab10Ambig_HiGini.pdf", height = 5, width = 5)
  plot(pca$perc_var)
  dev.off()
  
  df_out <- as.data.frame(matrix(nrow = length(rownames(forPCA)), ncol = 9))
  names(df_out) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "Name", "Extra_Name")
  df_out$PC1 <- pca$points[,1]
  df_out$PC2 <- pca$points[,2]
  df_out$PC3 <- pca$points[,3]
  df_out$PC4 <- pca$points[,4]
  df_out$PC5 <- pca$points[,5]
  df_out$PC6 <- pca$points[,6]
  df_out$PC7 <- pca$points[,7]
  df_out$Name <-  c(as.character(row.names(forPCA)))
  
  #This merges the file with the group information and the PC information
  df_out <- merge(df_out, GROUPS, by = "Name")
  
  df_out_controls <- subset(df_out, Data_Source == "Dawe_Lab_1" | Data_Source == "Dawe_Lab_2")
  df_out_controls$`RF Class` <- df_out_controls$Ab10_Status
  df_out_controls$Control <- "Yes"
  df_out_controls <- df_out_controls[,c("Name", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "Data_Source", "Accession", "Group", "Maize_Type", "Ab10_Status", "B_Chrom_Status", "Latitude", "Longitude", "Altitude", "RF Class", "Control")]
  
  
  df_out_exp <- subset(df_out, Data_Source != "Dawe_Lab_1" | Data_Source != "Dawe_Lab_2")
  df_out_exp <- merge(df_out_exp, df_col_Ab10Type, by = "Name")
  df_out_exp$Control <- "No"
  df_out_exp <- df_out_exp[,c("Name", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "Data_Source", "Accession", "Group", "Maize_Type", "Ab10_Status", "B_Chrom_Status", "Latitude", "Longitude", "Altitude", "RF Class", "Control")]
  
  df_out_edit <- rbind(df_out_controls, df_out_exp)
  
  Ab10_colors <- c("blue", "darkgreen", "magenta", "darkgrey")
  data_colors <- c("gold", "limegreen", "mediumorchid", "#008080", "#310062")
  shapes <- c(8, 15, 16, 3, 4, 7, 17, 18)
  
  #This loop plots the first two PCAs against each other
  pdf("PCA1_v_PCA2_Ab10Ambig_HiGini.pdf", height = 8, width = 15)
  ggplot() +
    geom_point(data=df_out_edit, alpha = 0.5, aes(x=PC1, y=PC2, color=`RF Class`, shape = Maize_Type, size = Control)) +
    geom_text(data=subset(df_out_edit, Ab10_Status == "Ab10-I"), aes(x=PC1, y=PC2, label=Name)) +
    ggtitle("PCA of GBS Tags\n on the Ab10") +
    scale_shape_manual(values= shapes) +
    scale_color_manual(values = Ab10_colors) +
    labs(color="Ab10", shape = "Maize Type") +
    labs(x="PC1 32%", y="PC2 18%") +
    theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
  dev.off()
  
  pdf("PCA1_v_PCA3_Ab10Ambig_HiGini.pdf", height = 8, width = 15)
  ggplot() +
    geom_point(data=df_out_edit, alpha = 0.5, aes(x=PC1, y=PC3, color=`RF Class`, shape = Maize_Type, size = Control)) +
    geom_text(data=subset(df_out_edit, Ab10_Status == "Ab10-I"), aes(x=PC1, y=PC3, label=Name)) +
    ggtitle("PCA of GBS Tags\n on the Ab10") +
    scale_shape_manual(values= shapes) +
    scale_color_manual(values = Ab10_colors) +
    labs(color="Ab10", shape = "Maize Type") +
    labs(x="PC1 32%", y="PC3 15%") +
    theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
  dev.off()
  
  pdf("PCA1_v_PCA4_Ab10Ambig_HiGini.pdf", height = 8, width = 15)
  ggplot() +
    geom_point(data=df_out_edit, alpha = 0.5, aes(x=PC1, y=PC4, color=`RF Class`, shape = Maize_Type, size = Control)) +
    geom_text(data=subset(df_out_edit, Ab10_Status == "Ab10-I"), aes(x=PC1, y=PC4, label=Name)) +
    ggtitle("PCA of GBS Tags\n on the Ab10") +
    scale_shape_manual(values= shapes) +
    scale_color_manual(values = Ab10_colors) +
    labs(color="Ab10", shape = "Maize Type") +
    labs(x="PC1 32%", y="PC4 11%") +
    theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
  dev.off()
  
  pdf("PCA1_v_PCA2_Ab10Ambig_HiGini_data.pdf", height = 8, width = 15)
  ggplot() +
    geom_point(data=df_out_edit, alpha = 0.5, aes(x=PC1, y=PC2, color=Data_Source, shape = Maize_Type, size = Control)) +
    geom_text(data=subset(df_out_edit, Ab10_Status == "Ab10-I"), aes(x=PC1, y=PC4, label=Name)) +
    ggtitle("PCA of GBS Tags\n on the Ab10") +
    scale_shape_manual(values= shapes) +
    scale_color_manual(values = data_colors) +
    labs(color="Ab10", shape = "Maize Type") +
    labs(x="PC1 32%", y="PC2 18%") +
    theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
  dev.off()
  
  pdf("PCA1_v_PCA2_Ab10Ambig_HiGini_data.pdf", height = 8, width = 15)
  ggplot() +
    geom_point(data=df_out_edit, alpha = 0.5, aes(x=PC1, y=PC2, color=Data_Source, shape = Maize_Type, size = Control)) +
    geom_text(data=subset(df_out_edit, Ab10_Status == "Ab10-I"), aes(x=PC1, y=PC4, label=Name)) +
    ggtitle("PCA of GBS Tags\n on the Ab10") +
    scale_shape_manual(values= shapes) +
    scale_color_manual(values = data_colors) +
    labs(color="Ab10", shape = "Maize Type") +
    labs(x="PC1 32%", y="PC2 18%") +
    theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
  dev.off()
  
  
  #This adds the Approximate Classes and plots out only the ambiguous samples. 
  df_out_edit2 <-  df_out_edit
  df_out_edit2 <- merge(df_out_edit, df_col_Ab10TypeAmbig, by="Name")
  
  pdf("PCA1_v_PCA2_Ab10Ambig_HiGini_Ambiguous.pdf", height = 8, width = 15)
  ggplot() +
    geom_point(data=df_out_edit2, aes(x=PC1, y=PC2, color=`Approximate Class`, shape = Maize_Type, alpha = Confidence, size = Confidence)) +
    geom_text(data=subset(df_out_edit2, Ab10_Status == "Ab10-I"), aes(x=PC1, y=PC2, label=Name)) +
    ggtitle("PCA of GBS Tags\n on the Ab10\n Ambiguous Samples Only") +
    scale_shape_manual(values= shapes) +
    scale_color_manual(values = Ab10_colors) +
    labs(color="Ab10", shape = "Maize Type") +
    labs(x="PC1 32%", y="PC2 18%") +
    theme(plot.title = element_text(hjust=0.5, size = 20), axis.title = element_text(size = 18), axis.text = element_text(size = 15), legend.title = element_text(size = 18), legend.text = element_text(size = 15), legend.position = "right")
  dev.off()
  
