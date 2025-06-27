#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)
file = args[1]
#module load R/4.4.1-foss-2022b

#Load the libraries
#install.packages("raster")
library(raster)

#Set the working directory with the worldclim2 data
setwd("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/ClimateData/SSP126")
dir <- "/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/ClimateData/SSP126"

#Drop the .tif extension
file <- gsub(".tif", "", file)

#Loop through all the .tif files and convert them to .asc files to be compatible with MaxEnt
print(paste("##### Working of file", file))
f <- paste(dir, "/", file, ".tif", sep="")
r <- raster(f)
ra <- aggregate(r, fact=2)  ## By default aggregates using mean, but see fun=
writeRaster(ra, paste(dir, "/", file, ".asc", sep=""), format="ascii")
#I manually moved all asc files to a folder labeled asc



