library(vroom)
library(data.table)

setwd("K10L2_Model/Experimental_K10L2N10Model/QualityControl")

#This cbinds all of the sorted columns together
IT=1
DATA <- vroom::vroom(paste("K10L2_QualityControl_KMeansGroupingControls_MixedControls", IT, "table", sep="."))


#This sorts all the files
for(i in 1:24) {
  a <- i+1
  TEMP <- vroom::vroom(paste("K10L2_QualityControl_KMeansGroupingControls_MixedControls", a, "table", sep="."))
  DATA <- rbind(DATA, TEMP)
}

#This writes out the final file
write.csv(DATA, file="K10L2_QualityControl_KMeansGroupingControls_MixedControls_All.csv", row.names = FALSE, quote = FALSE)
