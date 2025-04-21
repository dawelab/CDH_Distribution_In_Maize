library(vroom)
library(data.table)

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChrom_Model/Experimental_BChromModel_High/KmeansGroups_Files")

#This sorts all the files
for(i in 1:25) {
  DATA <- vroom::vroom(paste("Kmeans_Exp_Groups_High", i, "table", sep="."))
  DATA <- DATA[order(DATA$Name),]
  fwrite(DATA, file=paste("Kmeans_Exp_Groups_High_sorted", i, "table", sep="."))
}

#This cbinds all of the sorted columns together
IT=1
DATA <- vroom::vroom(paste("Kmeans_Exp_Groups_High_sorted", IT, "table", sep="."))

for(i in 1:24) {
  a=i+1
  TEMP <- vroom::vroom(paste("Kmeans_Exp_Groups_High_sorted", a, "table", sep="."))
  DATA <- cbind(DATA, TEMP[,-c(1)])
}

#This writes out the final file
write.csv(DATA, file="Kmeans_Exp_Groups_High_All.csv", row.names = FALSE, quote = FALSE)
