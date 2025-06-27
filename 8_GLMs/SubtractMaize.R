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
name <- paste0(type, "_projection_avg.suitable.minmax.extend")

#Identify the maize raster for this particular layer
if(model=="Current_Full") {
    model2 <- "SSPCurrent"
    dir2 <<- paste0("Maize_MaxEntStand4")
}

if(model=="Current_SSP") {
    model2 <- "SSPCurrent"
    dir2 <<- paste0("Maize_MaxEntStand4_", model2)
}

if(model!="Current_SSP" & model!="Current_Full") {
    dir2 <<- paste0("Maize_MaxEntStand4_", model)
}

#Load the Maize raster for this model
maize <- stack(paste0("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/", dir2, "/Maize_projection_avg.suitable.minmax.extend.asc"))

#load the CDH raster
cdh <- stack(paste0(name, ".asc"))

#This combines both rasters into a single stack
rast_stack <- stack(cdh,maize)
#This subtracts the maize raster from the CDH raster
options(scipen=999)
sub <- calc(rast_stack, function(x) { (x[1]-x[2])})

#df <- as.data.frame(sub, xy=TRUE)
#df_sub <- gsub("0", "NA", df)
#raster <- raster(df_sub)

#write the raster
writeRaster(sub, paste0(name, ".minusmaize.asc"), bylayer=TRUE, format="ascii", overwrite=TRUE)
