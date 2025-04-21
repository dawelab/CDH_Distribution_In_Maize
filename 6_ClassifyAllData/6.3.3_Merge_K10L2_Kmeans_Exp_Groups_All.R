library(vroom)
library(data.table)

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_K10L2Contigs/K10L2_Model/Experimental_K10L2N10Model/KmeansGroups_Files")

#This sorts all the files
for(i in 1:25) {
  DATA <- vroom::vroom(paste("K10L2_Kmeans_Exp_Groups_MixedControls", i, "table", sep="."))
  DATA <- DATA[order(DATA$Name),]
  fwrite(DATA, file=paste("K10L2_Kmeans_Exp_Groups_MixedControls_sorted", i, "table", sep="."))
}

#This cbinds all of the sorted columns together
IT=1
DATA <- vroom::vroom(paste("K10L2_Kmeans_Exp_Groups_MixedControls_sorted", IT, "table", sep="."))

for(i in 1:24) {
  a=i+1
  TEMP <- vroom::vroom(paste("K10L2_Kmeans_Exp_Groups_MixedControls_sorted", a, "table", sep="."))
  DATA <- cbind(DATA, TEMP[,-c(1)])
}

#This writes out the final file
write.csv(DATA, file="K10L2_Kmeans_Exp_Groups_MixedControls_All.csv", row.names = FALSE, quote = FALSE)

