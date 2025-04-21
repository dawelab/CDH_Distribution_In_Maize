library(vcfR)
library(stringr)
library(statgenGWAS)
install.packages("kinship2")
library(kinship2)

#Load the VCF File
VCF <- read.vcfR("/Users/user/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/CoordinateOnly_Ab10BChrom.chr1to9.filt4.1Perc.vcf.gz")

#Load the data with all of the CDH calls
GROUPS <- read.csv("~/University_of_Georgia/Dawe_Lab_Documents/Ab10_Global_Survey/R_Sessions/SGE_Distribution_Paper/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumCompleteWholeGenomePC.csv")

#This extracts only the genotypes
print("#####Extracting genotypes")
GT<- extract.gt(VCF)

#This converts the genotypes to an additive coding so that it can be used in a PCA
print("#####Converting genotypes to additive coding")
convert <- function(x) {
  y <- str_split((x), ":")[[1]][1]
  y <- sub("0/1", 1, y)
  y <- sub("1/0", 1, y)
  y <- sub("1/1", 2, y)
  y <- sub("0/0", 0, y)
  x <<-y
}
GT <- apply(GT, c(1, 2), convert)

#This creates a data frame with all the SNP positions
MAP <- data.frame(chr=c(rep(NA, length(rownames(GT)))), pos=rep(NA, length(rownames(GT))))
rownames(MAP) <- rownames(GT)

reformat_chr <- function(x) {
  chr <- str_split((x), "_")[[1]][1]
  chr <- gsub("S", "", chr)
  x <<-chr
}

reformat_pos <- function(x) {
  pos<- str_split((x), "_")[[1]][2]
  x <<- pos
}

MAP$chr <- lapply(rownames(MAP), reformat_chr)
MAP$chr <- as.numeric(MAP$chr)
MAP$pos <- lapply(rownames(MAP), reformat_pos)
MAP$pos <- as.numeric(MAP$pos)

#This makes the response variable (ie phenotype)
CDH <-   GROUPS[,c("Name", "KMeans_Ab10", "KMeans_BChrom", "KMeans_K10L2")]
colnames(CDH) <- c("genotype", "KMeans_Ab10", "KMeans_BChrom", "KMeans_K10L2")

#This reformats all the predictors as a numerical
CDH$KMeans_Ab10 <- gsub("Ab10", 1, CDH$KMeans_Ab10)
CDH$KMeans_Ab10 <- gsub("N10", 0, CDH$KMeans_Ab10)
CDH$KMeans_Ab10 <- gsub("Ambiguous", NA, CDH$KMeans_Ab10)

CDH$KMeans_BChrom <- gsub("Yes", 1, CDH$KMeans_BChrom)
CDH$KMeans_BChrom <- gsub("No", 0, CDH$KMeans_BChrom)
CDH$KMeans_BChrom <- gsub("Ambigious", NA, CDH$KMeans_BChrom)

CDH$KMeans_K10L2 <- gsub("K10L2", 1, CDH$KMeans_K10L2)
CDH$KMeans_K10L2 <- gsub("N10", 0, CDH$KMeans_K10L2)
CDH$KMeans_K10L2 <- gsub("Ambigious", NA, CDH$KMeans_K10L2)

CDH$KMeans_Ab10 <- as.numeric(CDH$KMeans_Ab10)
CDH$KMeans_BChrom <- as.numeric(CDH$KMeans_BChrom)
CDH$KMeans_K10L2 <- as.numeric(CDH$KMeans_K10L2)

#This generates the full data object 
#This converts the matrix to a dataframe
GT <-t(GT)
OBJ <- createGData(geno = GT, map = MAP, pheno = CDH)
RECODE <- codeMarkers(OBJ, impute = TRUE)

print(RECODE$markers)

KIN <- statgenGWAS::kinship(RECODE)
print(OBJ)
