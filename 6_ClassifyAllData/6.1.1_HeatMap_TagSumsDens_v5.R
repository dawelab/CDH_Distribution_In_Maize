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
setwd("/Volumes/Transcend")

#setwd("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_v7")

########################################This loads and preps the data 
GROUPS <- vroom::vroom("Swarts_AllControls_RomeroNavarro_Groups_Env.csv")

MERGE_Ab10Hap_RPM <- vroom::vroom("BWAaln_SwartsAllControlsRN_v_Ab10HIFIBChrom.Ab10Hap.TagEdit.TagIDLoc.RPM.table")

#28 columns were classed as logicals because they have all NAs, these are the ones that had 0 coverage

#This drops any column that had an average coverage of 0, resulting in all values being NA
DT <- as.data.table(MERGE_Ab10Hap_RPM)
DT <- DT[,which(unlist(lapply(DT, function(x)!all(is.na(x))))),with=F]
MERGE_Ab10Hap_RPM <- as.data.frame(DT)

#This loads in the BINS file
BINS <- read_excel("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/Bins_NoOverlap.table.xlsx")

#This removes any tag with a MAPQ less than 20
MERGE_Ab10Hap_RPM_FILT <- subset(MERGE_Ab10Hap_RPM, BWA_MAPQ >= 20)

#this adds a bin value to each line in the MERGE_Ab10Hap_RPM file 
MERGE_Ab10Hap_RPM_FILT$bin <- NA
i=1
x=1
for (i in 1:nrow(MERGE_Ab10Hap_RPM_FILT)) {
  print(i)
  for (x in 1:nrow(BINS)) {
    BIN_NUM <- BINS[x,4][[1]]
    START <- BINS[x,2][[1]]
    END <- BINS[x,3][[1]]
    ####CHANGE THIS NUMBER TO 3 WHEN YOU FIX THE PROBLEM WITH THE START
    if (MERGE_Ab10Hap_RPM_FILT[i,3] > START & MERGE_Ab10Hap_RPM_FILT[i,3] < END) {
      MERGE_Ab10Hap_RPM_FILT[i,ncol(MERGE_Ab10Hap_RPM_FILT)] <- BIN_NUM
    }
  }
}

#This outputs the bins as well as the tag coordinates to verify correct identification
test <- MERGE_Ab10Hap_RPM_FILT[,c(1:6,ncol(MERGE_Ab10Hap_RPM_FILT))]

#This generates a data frame with of all the line names found in the data sest
Name <- as.data.frame(colnames(MERGE_Ab10Hap_RPM_FILT[-c(1:7)]))
colnames(Name) <- "Name"

#This section alters the names of the Romero-Navarro data so that it can be merged with the groups file
Name$Extra_Name <- NA
Name$OG_Name<- Name$Name
for(i in 1:nrow(Name)){
  NAME <- Name[i,1]
  NAME <- as.character(NAME)
  if(startsWith(NAME, "SEEDGWAS") == TRUE){
    TEMP <- str_split(NAME, ":", n=2)
    Name[i,1] <- TEMP[[1]][1]
    Name[i,2] <- TEMP[[1]][2]
  }
  if(startsWith(NAME, "12S00") == TRUE){
    TEMP <- str_split(NAME, ":", n=2)
    Name[i,1] <- TEMP[[1]][1]
    Name[i,2] <- TEMP[[1]][2]
  }
  if(startsWith(NAME, "RIMMA") == TRUE){
    TEMP <- str_split(NAME, "_", n=2)
    Name[i,1] <- TEMP[[1]][1]
    Name[i,2] <- TEMP[[1]][2]
  }
}

#This section isolates the Romero Navarro Data and merges it with the groups.
Name_RM <- subset(Name, is.na(Extra_Name) == FALSE)
GROUPS_RM <- subset(GROUPS, Data_Source == "Romero-Navarro_etal_2017")
Name_RM <- merge(Name_RM, GROUPS_RM, by = "Name")
Name_RM$CleanName <- paste(Name_RM$Name, ".RM", sep="")

#This section merges everything else that is not Romero-Navarro Data
Name_ELSE <- subset(Name, is.na(Extra_Name) == TRUE)
GROUPS_ELSE <- subset(GROUPS, Data_Source != "Romero-Navarro_etal_2017")
Name_ELSE <- merge(Name_ELSE, GROUPS_ELSE, by = "Name")
#You really do want this to be the original name
Name_ELSE$CleanName <- Name_ELSE$Name

#This section specifically extracts the RIMMA and blank data from the swarts data sets, which is lost in the other two filters
Name_RA <- subset(Name, startsWith(Name, "RIMMA") == TRUE | Name == "BLANK")
Name_RA <- merge(Name_RA, GROUPS, by = "Name")
Name_RA$CleanName <- Name_RA$OG_Name

#This merges all of the data together
df_extra <- rbind(Name_RM, Name_ELSE, Name_RA)

#This pulls out everything that is not Romero-Navarro
df_extra_OTH <- subset(df_extra, Data_Source != "Romero-Navarro_etal_2017")

##################This section isolates all of the data that does not need to be altered
#This identifies all the index numbers for non Romero-Navarro Lines
col.num <- which(colnames(MERGE_Ab10Hap_RPM_FILT) %in% df_extra_OTH$OG_Name)
#This appends the meta data columns to this list
val <- c(1:6,18224)
col.num <- append(val, col.num)
#This selects all of these columns
MERGE_Ab10Hap_RPM_FILT_OTH <- MERGE_Ab10Hap_RPM_FILT[,c(col.num)]

##################Sums the columns for RM
#This creates a list of all the unique names from the Romero-Navarro dataset
SAMPS <- unique(Name_RM$Name)

###################This section takes the mean of all the rows for all the Romero-Navarro Data originating from the same individual

#This creates a list of all the unique names from the Romero-Navarro dataset. This was already run above
SAMPS <- unique(Name_RM$Name)

i=1
for(i in 1:length(SAMPS)) {
  x <- SAMPS[i]
  index <- which(Name_RM$Name %in% x)
  OGNAME <- Name_RM[index, "OG_Name"]
  index2 <- which(colnames(MERGE_Ab10Hap_RPM_FILT) %in% OGNAME)
  TEMP1 <- as.data.frame(MERGE_Ab10Hap_RPM_FILT[,index2], drop = FALSE)
  MEANS <- as.data.frame(rowMeans(TEMP1))
  CLEAN_NAME <- Name_RM[index, "Name"][1]
  colnames(MEANS) <- paste(CLEAN_NAME,".RM", sep="")
  MERGE_Ab10Hap_RPM_FILT_OTH <<- cbind(MERGE_Ab10Hap_RPM_FILT_OTH, MEANS)
}

MERGE_Ab10Hap_RPM_FILT_FIX <- MERGE_Ab10Hap_RPM_FILT_OTH

fwrite(MERGE_Ab10Hap_RPM_FILT_FIX, file="BWAaln_SwartsAllControlsRN_v_Ab10HIFIBChrom.Ab10Hap.TagEdit.TagIDLoc.RPM.RMMean.table", row.names = FALSE)

#This evaluates how much missing data is in each sample
MERGE_Ab10Hap_RPM_FILT_FIX_SHARE <- subset(MERGE_Ab10Hap_RPM_FILT_FIX, (start >= 152050000 & end <= 156350000) | (start >= 158250000 & end <= 166820000))

MISS <- data.frame(Name = colnames(MERGE_Ab10Hap_RPM_FILT_FIX_SHARE[-c(1:7)]), Missing = colSums(MERGE_Ab10Hap_RPM_FILT_FIX_SHARE[-c(1:7)] == 0))

MISS$Perc_Missing <- (MISS$Missing/nrow(MERGE_Ab10Hap_RPM_FILT_FIX_SHARE))*100

#This plots the percentage of missing tags per sample
ggplot(data=MISS, aes(x=Perc_Missing)) +
  geom_histogram(bins = 100) +
  scale_x_continuous(breaks=c(1:100)) +
  ggtitle("Percent Missing Data in Shared Region")

ggsave("PercentSharedRegionMissingData_Histogram.png")

#This drops any sample with more than 99% missing data, I may want to change this filter later
MISS_DROP <- subset(MISS, Perc_Missing >= 99)
col.num <- which(colnames(MERGE_Ab10Hap_RPM_FILT_FIX) %in% MISS_DROP$Name)
MERGE_Ab10Hap_RPM_FILT_FIX_D <- MERGE_Ab10Hap_RPM_FILT_FIX[,-c(col.num)]

#This function goes over each row and divides each value by the max in that row 
MinMax = function(xx) { sweep(xx, 1, apply(xx, 1, max), '/') }

######## Tag Sums
#This sums the number of tag numbers across bins, but excludes the bin column from the summing

A_SUM = as.data.frame(apply(MERGE_Ab10Hap_RPM_FILT_FIX_D[,-c(1:7)], 2, function(xx) { by(xx, MERGE_Ab10Hap_RPM_FILT_FIX_D$bin, sum)  }))

fwrite(A_SUM, file ="TagSums_AllSamples.csv")

######## Tag Density
#This sums the number of tag numbers across bins, but excludes the bin column from the summing
A_DENS = as.data.frame(apply(MERGE_Ab10Hap_RPM_FILT_FIX_D[,-c(1:7)], 2, function(xx) { by(xx, MERGE_Ab10Hap_RPM_FILT_FIX_D$bin, function(yy) { sum(yy > 0) }) }))

fwrite(A_SUM, file ="TagDens_AllSamples.csv")

#This takes the square root of each value

B_DENS <- sqrt(A_DENS)

#This combines the density and sum data
B_COMB <- A_SUM+B_DENS

fwrite(B_COMB, file ="TagSum_Plus_TagDens_AllSamples.csv")

#This scales the data so that the min is always 0 and the max is always 1. 
C <- MinMax(B_COMB)

#This adds bins back
C$bin <- rownames(C)

C_bed <- merge(BINS, C, by="bin")

#This drops the identifier columns to create a numeric matrix
rownames(C_bed) <- C_bed$bin
heat_data <- C_bed[,-c(1:4)]

#This makes a data frame for column annotation

#This seleects only unique entries, columns 2 and 3 include the extra Romero Navarro names to differentiate each replicate and need to be removed to find unique lines
df_extra_u <- unique(df_extra[,-c(2:3)])

#Select only the relevant columns
df_col <- df_extra_u[,c("CleanName", "Ab10_Status", "Data_Source")]

#This removes the samples dropped above for low coverage
row.num <- which(df_col$CleanName %in% MISS_DROP$Name)
#This selects all of these columns
df_col <- df_col[-c(row.num),]
#There is an extra BLANK, I'm not sure where it comes from, but it needs to be removed 
df_col <- subset(df_col, CleanName != "BLANK")

####This is to make it so the RIMMA names are are compatible in both datasets
#i=1
#for(i in 1:nrow(df_col)) {
#  df_col[i,1] <- strsplit(df_col[i,1],"[.]")[[1]][1]
#}

rownames(df_col) = df_col$CleanName
df_col<- df_col[,2:3, drop=FALSE]
colnames(df_col) <- c("Chr10 Type", "Data Source")
df_col$`Chr10 Type` <- as.factor(df_col$`Chr10 Type`)
df_col$`Data Source` <- as.factor(df_col$`Data Source`)

#This section reorders the data to match the heat data
df_col <- df_col[order(df_col$`Data Source`),]

#This line reorders the heat_data to match the df_col
heat_data <- heat_data[,rownames(df_col)]

#This line sets the rownames in the df as the column names in the heat map data. IT is very important
#rownames(df_col[1:3,]) = colnames(heat_data2[,1:3])

#This creates a dataframe for the position data 
pos <- C_bed$start
df_row <- as.data.frame(pos)
df_row$feature <- " "
i=1
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
summary(df_col$`Chr10 Type`)

colors <- list(`Chr10 Type` = c("Ab10-I" = "black","Ab10-II" = "grey20","Ab10-III" =  "grey40", "Ab10_Unknown" = "navy", "K10L2" = "grey60", "N10" ="grey80", "Unknown" = "white"), Feature = c('TR1' = "#56B4E9", 'trkin' = "#0072B2", 'Shared region' = "#E69F00", 'knob 180' = "#D55E00", 'kindr' = "#CC79A7", 'kin10-like' =  "#009E73", " " = "white"), `Data Source` = c("Dawe_Lab_1" = "gold", "Dawe_Lab_2" = "chartreuse2", "Romero-Navarro_etal_2017" = "#008080", "Swarts_etal_2017" = "#310062", "NA" = "white"))

pdf("HeatMap_TagSumsSqrtTagDens_MinMax_All.pdf", height = 10, width = 10)
pheatmap(heat_data, scale = "none", show_rownames = TRUE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE, treeheight_row = 0, treeheight_col = 0, annotation_row = df_row, annotation_colors= colors, color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrBr")))(500)), annotation_col = df_col, border_color = NA, main="Chr10 Haplotypes\nTag Index 1MB bins")
dev.off()

pdf("HeatMap_TagSumsSqrtTagDens_MinMax_ColClust_All.pdf", height = 10, width = 10)
pheatmap(heat_data, scale = "none", show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = TRUE, clustering_method = "ward.D", treeheight_row = 0, treeheight_col = 50, annotation_row = df_row, annotation_colors= colors, color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrBr")))(100)), annotation_col = df_col, border_color = NA, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_TagSumsSqrtTagDens_MinMax_ColClust_RowClust_All.pdf", height = 10, width = 20)
pheatmap(heat_data, scale = "none", show_rownames = FALSE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D", treeheight_row = 50, treeheight_col = 50, annotation_row = df_row, annotation_colors= colors, color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrBr")))(100)), annotation_col = df_col, border_color = NA, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()

pdf("HeatMap_TagSumsSqrtTagDens_MinMax_ColClust_RowClust_RowNames_All.pdf", height = 10, width = 10)
pheatmap(heat_data, scale = "none", show_rownames = TRUE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D", treeheight_row = 50, treeheight_col = 50, annotation_row = df_row, annotation_colors= colors, color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrBr")))(100)), annotation_col = df_col, fontsize_row = 3, border_color = NA, main="Chr10 Haplotypes\nTag Index in 1MB bins")
dev.off()

res <- pheatmap(heat_data, scale = "none", show_rownames = TRUE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D", treeheight_row = 50, treeheight_col = 50, annotation_colors= colors, color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrBr")))(100)), annotation_col = df_col, fontsize_row = 3, border_color = NA, main="Chr10 Haplotypes\nTag Index in 1MB bins")

Ab10_clusts <- as.data.frame(cbind(colnames(heat_data), cluster = cutree(res$tree_col,k = 2)))
colnames(Ab10_clusts) <- c("Name", "Ab10 Assignment")
Ab10_Yes <- subset(Ab10_clusts, Ab10_clusts$`Ab10 Assignment` == 1)

Ab10_clusts$`Ab10 Assignment` <- gsub(1, "Yes", Ab10_clusts$`Ab10 Assignment`)
Ab10_clusts$`Ab10 Assignment` <- gsub(2, "No", Ab10_clusts$`Ab10 Assignment`)

fwrite(Ab10_clusts, file ="Ab10_Assignments_AllSamples_v1.csv")
