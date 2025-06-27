
#module load R/4.4.1-foss-2022b
library(raster)
library(sf)

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/MaxEnt_EnvData_Subset/SSP370_projection")
FILES <- c("wc2.1_30s_elev.stand.asc")

#This loads in one world clim2 layer
TEMPLATE <- stack("wc2.1_30s_bio_18.stand.asc")

#This projects all of the other layers to the same projection, resolution, and extent as the template
for(i in FILES) {
  stack <- stack(i)
  #This diasables printing in scientific notation. For some reason it is required
  options(scipen=999)
  FIX <- projectRaster(stack, TEMPLATE)
  writeRaster(FIX, paste0(names(FIX),"2.asc"), bylayer=TRUE, format="ascii", overwrite=TRUE)
}