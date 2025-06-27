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
name <- paste0(type, "_projection_avg.suitable.minmax")

#Load the raster
r <- stack(paste0(name,".asc"))

#This changes the extent to the extent of the original environemntal layers
TEMPLATE <- stack("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/ClimateData/Current/asc/wc2.1_30s_bio_18.asc")

#This diasables printing in scientific notation. For some reason it is required
options(scipen=999)
FIX <- projectRaster(r, TEMPLATE)

#This changes all the newly added NA values to 0 so that I can combine each CDH raster with the maize raster
# Fill all missing values with 0
FIX[is.na(FIX)] <- 0 

#This writes the altered raster to a new file
writeRaster(FIX, paste0(type, "_projection_avg.suitable.minmax.extend.asc"), bylayer=TRUE, format="ascii", overwrite=TRUE)