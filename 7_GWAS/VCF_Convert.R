library(stringr)
library(vroom)
library(data.table)

setwd("/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_v6")

VCF <- vroom::vroom("AllData_v1_Ab10Spec.vcf.recode.format.vcf")

convert <- function(x) {
  y <- str_split((x), ":")[[1]][1]
  y <- sub("0/1", "2", y)
  y <- sub("1/1", "2", y)
  y <- sub("0/0", "1", y)
  if(grepl(".", y, fixed = TRUE)==TRUE) {
    y<- 0
  }
  x <<-y
}

VCF[,-c(1:9)] <- apply(VCF[,-c(1:9)], c(1, 2), convert)

fwrite(VCF, file = "AllData_v1_Ab10Spec.vcf.recode.format.convert.vcf")
