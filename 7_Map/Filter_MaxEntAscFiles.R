#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)
dir = args[1]
type = args[2]
model = args[3]

#dir <- "Ab10_MaxEntStand4"
#type <- "Ab10"
#model <- "Current_Full"

#module load R/4.4.1-foss-2022b
library(raster)
library(sp)

#Set the working directory 
setwd(paste0("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/",dir))

#Extract the suitability threshold value from the maxent csv file
csv <- read.csv("maxentResults.csv")
threshold <- csv$Maximum.test.sensitivity.plus.specificity.Logistic.threshold[11]

#This prepares the correct name format for the model
if(model != "Current_Full" & model != "Current_SSP") {
    name <<- paste0(type, "_", model, "_projection_avg.asc")
} 
if(model == "Current_Full" | model == "Current_SSP") {
    name <<- paste0(type, "_projection_avg.asc")
} 

#Load the raster
r <- stack(name)

#This filters the raster to only include values greater than the threshold
r[r[] < threshold ] = NA

#This calculates the total area of the suitable habitat ignoring NA values
area <- area(r, na.rm=TRUE)

#This sums all the areas
Sum_area <- cellStats(area, 'sum')

#This reads in the summary table
table <- read.csv("/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Stand4_Area_Table.csv")

#This identifies the correct column and row in the summary table
row <- which(table$Model == model)
col <- which(colnames(table) == type)

#This assigns the square meter value to the correct cell in the summary table
table[row, col] <- Sum_area

#This writes out the updated summary table
write.csv(table, "/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Stand4_Area_Table.csv", row.names=FALSE)

#This writes the filtered raster to a new file
writeRaster(r, paste0(type, "_projection_avg.suitable.asc"), bylayer=TRUE, format="ascii", overwrite=TRUE)
