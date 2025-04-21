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
name <<- paste0(type, "_projection_avg.suitable")

#Load the raster
r <- stack(paste0(name,".asc"))

#This identifies the max value in the raster
mean <- cellStats(r, 'mean')
#This identifies the min value in the raster
sd <- cellStats(r, 'sd')
#This min/max scales the raster
calc(r, fun = function(x) ((x-mean[1][[1]])/sd[1][[1]]), filename=paste0(name, ".norm.asc"), overwrite=TRUE)

#This loads the normalized raster
r <- stack(paste0(name, ".norm.asc"))

#This checks that the max value is 1 and the min is 0
max <- cellStats(r, 'max')
print("This is the max:")
print(max)
min <- cellStats(r, 'min')
print("This is the min:")
print(min)