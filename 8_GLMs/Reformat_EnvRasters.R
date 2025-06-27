#module load R/4.4.1-foss-2022b
library(raster)
library(sf)

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/ClimateData/Current/asc")
FILES <- c("wc2.1_30s_vapr_avg.asc")

#this reads in the shapefile with the area of interest
AOI <- read_sf("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/MaxEnt_EnvData_Subset/AOI_ConvexHull_CountryBoundaries.shp")

#This crops one of the world clim layers to only the americas
i="wc2.1_30s_bio_18.asc"
stack <- stack(i)
Cropped <- crop(x = stack, y = AOI)
writeRaster(Cropped, paste0("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/MaxEnt_EnvData_Subset/stand/", names(Cropped),".stand.asc"), bylayer=TRUE, format="ascii", overwrite=TRUE)

#This loads in the cropped worldclim2 layer as a template
TEMPLATE <- stack("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/MaxEnt_EnvData_Subset/stand/wc2.1_30s_bio_18.stand_1.asc")

#This projects all of the other layers to the same projection, resolution, and extent as the template
for(i in FILES) {
  stack <- stack(i)
   #This diasables printing in scientific notation. For some reason it is required
  options(scipen=999)
  FIX <- projectRaster(stack, TEMPLATE)
  writeRaster(FIX, paste0("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/MaxEnt_EnvData_Subset/stand/",names(FIX),".stand.asc"), bylayer=TRUE, format="ascii", overwrite=TRUE)
}


#setwd("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/SoilData")
#FILES <- c("sq5.asc")

#This projects all of the other layers to the same projection, resolution, and extent as the template
#for(i in FILES) {
#  stack <- stack(i)
#   #This diasables printing in scientific notation. For some reason it is required
#  options(scipen=999)
#  FIX <- projectRaster(stack, TEMPLATE)
#  writeRaster(FIX, paste0("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/MaxEnt_EnvData_Subset/stand/",names(FIX),".stand.asc"), bylayer=TRUE, format="ascii", overwrite=TRUE)
#}
