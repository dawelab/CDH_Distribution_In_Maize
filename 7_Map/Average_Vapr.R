#!/usr/bin/Rscript

#module load R/4.4.1-foss-2022b
library(raster)
library(sf)

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/ClimateData/Current/asc")
FILES <- c("wc2.1_30s_vapr_01.asc", "wc2.1_30s_vapr_02.asc", "wc2.1_30s_vapr_03.asc", "wc2.1_30s_vapr_04.asc", "wc2.1_30s_vapr_05.asc","wc2.1_30s_vapr_06.asc", "wc2.1_30s_vapr_07.asc", "wc2.1_30s_vapr_08.asc", "wc2.1_30s_vapr_09.asc", "wc2.1_30s_vapr_10.asc", "wc2.1_30s_vapr_11.asc", "wc2.1_30s_vapr_12.asc")

#This reads in the file
STACK <- stack(FILES)
#This averages the layers
AVG <- calc(STACK, fun = mean)
#This writes out the raster
writeRaster(AVG, "wc2.1_30s_vapr_avg.asc", bylayer=TRUE, format="ascii", overwrite=TRUE)
