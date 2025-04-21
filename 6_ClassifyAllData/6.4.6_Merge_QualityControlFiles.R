library(vroom)
library(data.table)

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChrom_Model/Experimental_BChromModel_Low/QualityControl")

#This cbinds all of the sorted columns together
IT=1
DATA <- vroom::vroom(paste("QualityControl_KMeansGroupingControls_Low", IT, "table", sep="."))


#This sorts all the files
for(i in 1:24) {
  a <- i+1
  TEMP <- vroom::vroom(paste("QualityControl_KMeansGroupingControls_Low", a, "table", sep="."))
  DATA <- rbind(DATA, TEMP)
}

#This writes out the final file
write.csv(DATA, file="QualityControl_KMeansGroupingControls_Low_All.csv", row.names = FALSE, quote = FALSE)
