#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)
dir = args[1]
file = args[2]
#module load R/4.4.1-foss-2022b

#Load the libraries
library(raster)

#Set the working directory with the worldclim2 data
setwd(dir)
#Load in the file as a stack
s1 <- stack(file)
#write out each layer as a .tif file
writeRaster(s1, paste0(names(s1),".tif"), bylayer=TRUE, format="GTiff", overwrite=TRUE)