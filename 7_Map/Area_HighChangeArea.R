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

#load the CDH raster
r <- stack(paste0(name, ".asc"))

#This filters the raster to only include values that are not 0
r[r[] == 0 ] = NA

#Calculate the mean and standard deviation
mean <- cellStats(r, "mean", na.rm=TRUE)
print("This is the mean:")
print(mean)

sd <- cellStats(r, "sd", na.rm=TRUE)
print("This is the sd:")
print(sd)

#This reports the raster max and min values
#max <- maxValue(r)
#print("This is the max value:")
#print(max)

#min <- minValue(r)
#print("This is the min value:")
#print(min)

#This calculates the values 2 standard deviations above the mean and below the mean
above <- mean + (2*sd)
below <- mean - (2*sd)

#This renmaes the raster
#tails <- r

#Filter the raster to only include values that are 2 standard deviations or more outside the mean
#tails[tails[] > below & tails[] < above ] = NA

#This writes the filtered raster to a new file
#writeRaster(tails, paste0(type, "_projection_avg.suitable.minmax.extend.minusmaize.2sd.asc"), bylayer=TRUE, format="ascii", overwrite=TRUE)

#This reads in the filtered raster
tails <- stack(paste0(type, "_projection_avg.suitable.minmax.extend.minusmaize.2sd.asc"))

# #This renmaes the raster
# lefttail <- tails

# #This isolates the below 2sd raster and calculates area
# lefttail[lefttail[] > below ] = NA

# #This writes out the below raster
# writeRaster(lefttail, paste0(type, "_projection_avg.suitable.minmax.extend.minusmaize.below2sd.asc"), bylayer=TRUE, format="ascii", overwrite=TRUE)

# #This calculates the area the the CDH is less likely to occur than maize
# below_area <- area(lefttail, na.rm=TRUE)
# below_area_sum <- cellStats(below_area, 'sum')

# #This reads in the summary table 
# below_table <- read.csv("/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Stand4_Below2SDMaize_Area_Table.csv")

# #This identifies the correct column and row in the summary table. The two tables are laid out in the same way so they can share a value
# row <- which(below_table$Model == model)
# col <- which(colnames(below_table) == type)

# #This assigns the square meter value to the correct cell in the summary table
# below_table[row, col] <- below_area_sum

# #This writes out the updated summary table
# write.csv(below_table, "/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Stand4_Below2SDMaize_Area_Table.csv", row.names=FALSE)

#This renmaes the raster
righttail <- tails

#This isolates the above 2sd raster and calculates area
righttail[righttail[] < above ] = NA

#This writes out the below raster
writeRaster(righttail, paste0(type, "_projection_avg.suitable.minmax.extend.minusmaize.above2sd.asc"), bylayer=TRUE, format="ascii", overwrite=TRUE)

#This calculates the area the the CDH is more likely to occur than maize
above_area <- area(righttail, na.rm=TRUE)
above_area_sum <- cellStats(above_area, 'sum')

#This reads in the summary table 
above_table <- read.csv("/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Stand4_Above2SDMaize_Area_Table.csv")

#This identifies the correct column and row in the summary table. The two tables are laid out in the same way so they can share a value
row <- which(below_table$Model == model)
col <- which(colnames(below_table) == type)

#This assigns the square meter value to the correct cell in the summary table
above_table[row, col] <- above_area_sum

#This writes out the updated summary table
write.csv(above_table, "/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Stand4_Above2SDMaize_Area_Table.csv", row.names=FALSE)
