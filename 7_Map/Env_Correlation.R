#module load R/4.4.1-foss-2022b
#This loads the libraries
library("raster")

#This loads the GPS data for the Swarts Data Set
SAMPS <- read.csv("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv", header=TRUE)
row.names(SAMPS) <- SAMPS$name
SAMPS_GPS <- SAMPS[,c("Longitude", "Latitude")] 
SAMPS_GPS$Longitude <- as.numeric(SAMPS_GPS$Longitude)
SAMPS_GPS$Latitude <- as.numeric(SAMPS_GPS$Latitude)


#This sets the working directory
setwd("/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/ClimateData/Current/asc")

#This lists all of the worldclim2 variabes to import
WC <- c("tmin", "tmax", "tavg", "prec", "srad", "wind", "vapr", "bio", "elev")

#This for loop reads in all of the worldclim2 variables as rastar files 
for(VAR in WC) {
  NAME <- (gsub(" ", "", paste(VAR, ".",  "ras")))
  list <- list.files(path="/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/ClimateData/Current/asc", pattern=VAR, full.names=TRUE)
  stack <- stack(list)
  assign(NAME, stack)
}

#This lists all of the worldclim2 variabes to import
SOIL <- c("sq1", "sq2", "sq3", "sq4", "sq5", "sq6", "sq7")

#This for loop reads in all of the worldclim2 variables as rastar files 
for(VAR in SOIL) {
  NAME <- (gsub(" ", "", paste(VAR, ".",  "ras")))
  list <- list.files(path="/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/SoilData", pattern=VAR, full.names=TRUE)
  stack <- stack(list)
  assign(NAME, stack)
}


#These commands extract the data for the GPS points represented in my sample files
tmin.data <- extract(tmin.ras, SAMPS_GPS)
tmax.data <- extract(tmax.ras, SAMPS_GPS)
tavg.data <- extract(tavg.ras, SAMPS_GPS)
prec.data <- extract(prec.ras, SAMPS_GPS)
srad.data <- extract(srad.ras, SAMPS_GPS)
wind.data <- extract(wind.ras, SAMPS_GPS)
vapr.data <- extract(vapr.ras, SAMPS_GPS)
bio.data <- extract(bio.ras, SAMPS_GPS)
elev.data <- extract(elev.ras, SAMPS_GPS)
nutavail.data <- extract(sq1.ras, SAMPS_GPS)
nutret.data <- extract(sq2.ras, SAMPS_GPS)
rootcond.data <- extract(sq3.ras, SAMPS_GPS)
oavail.data <- extract(sq4.ras, SAMPS_GPS)
salts.data <- extract(sq5.ras, SAMPS_GPS)
toxic.data <- extract(sq6.ras, SAMPS_GPS)
work.data <- extract(sq7.ras, SAMPS_GPS)

#This brings together all of my sample information with all of the WorldClim2 environmsntal data 
WC2.samps <- cbind(SAMPS, tmin.data, tmax.data, tavg.data, prec.data, srad.data, wind.data, vapr.data, bio.data, elev.data, nutavail.data, nutret.data, rootcond.data, oavail.data, salts.data, toxic.data, work.data)

#This writes the data set to a file
write.csv(WC2.samps, "Controls_Swarts_RomeroNavarro_Romay_Groups_Env.csv")

