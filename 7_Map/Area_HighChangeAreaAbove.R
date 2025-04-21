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

#This prepares the correct name format for the model
name <- paste0(type, "_projection_avg.suitable.minmax.extend.minusmaize")

#This loads the raster containing only values above 2 standard deviations of the mean
righttail <- raster(paste0(type, "_projection_avg.suitable.minmax.extend.minusmaize.above2sd.asc"))

#This calculates the area the the CDH is more likely to occur than maize
above_area <- area(righttail, na.rm=TRUE)
above_area_sum <- cellStats(above_area, 'sum')

#This reads in the summary table 
above_table <- read.csv("/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Stand4_Above2SDMaize_Area_Table.csv")

#This identifies the correct column and row in the summary table. The two tables are laid out in the same way so they can share a value
row <- which(above_table$Model == model)
col <- which(colnames(above_table) == type)

#This assigns the square meter value to the correct cell in the summary table
above_table[row, col] <- above_area_sum

#This writes out the updated summary table
write.csv(above_table, "/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Stand4_Above2SDMaize_Area_Table.csv", row.names=FALSE)
