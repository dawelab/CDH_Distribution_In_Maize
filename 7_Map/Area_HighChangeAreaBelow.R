#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)
dir = args[1]
type = args[2]
model = args[3]

#module load R/4.4.1-foss-2022b
library(raster)
library(sp)
library(spatialEco)

#Set the working directory 
setwd(paste0("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/",dir))

#This loads the raster containing only values below 2 standard deviations of the mean
lefttail <- raster(paste0(type, "_projection_avg.suitable.minmax.extend.minusmaize.below2sd.asc"))

#This calculates the area the the CDH is more likely to occur than maize
below_area <- area(lefttail, na.rm=TRUE)
below_area_sum <- cellStats(below_area, 'sum')

#This reads in the summary table 
below_table <- read.csv("/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Stand4_Below2SDMaize_Area_Table.csv")

#This identifies the correct column and row in the summary table. The two tables are laid out in the same way so they can share a value
row <- which(below_table$Model == model)
col <- which(colnames(below_table) == type)

#This assigns the square meter value to the correct cell in the summary table
below_table[row, col] <- below_area_sum

#This writes out the updated summary table
write.csv(below_table, "/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Stand4_Below2SDMaize_Area_Table.csv", row.names=FALSE)
