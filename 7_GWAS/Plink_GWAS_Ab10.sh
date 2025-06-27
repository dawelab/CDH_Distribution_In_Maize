#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Plink_GWAS_Ab10
#SBATCH --output Plink_GWAS_Ab10.out
#SBATCH --mem=50GB
#SBATCH --time=2:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load PLINK/1.9b_6.21-x86_64
module load BCFtools/1.15.1-GCC-11.3.0

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"
PHENO="/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/6_ClassifyAllData/6.5_Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumComplete.csv"

#This concatenates the imputed files
cd $DIR

bcftools concat -o $DIR/GWAS_Mo17_Ab10.chr.landrace.filt.imputed.vcf.gz $DIR/GWAS_Mo17_Ab10.chr.landrace.filt_chr*.imputed.vcf.gz

NAME="GWAS_Mo17_Ab10.chr.landrace.filt.imputed"

#This extracts the CDH call columns from the Pheno file
awk -F, '{print $1, $15, $17, $18, $20}' $PHENO > $DIR/CDHCall_ForPLINK.txt

#This changes the CDH coding
sed -i 's/Ab10/1/g' $DIR/CDHCall_ForPLINK.txt
sed -i 's/K10L2/1/g' $DIR/CDHCall_ForPLINK.txt
sed -i 's/Yes/1/g' $DIR/CDHCall_ForPLINK.txt
sed -i 's/No/0/g' $DIR/CDHCall_ForPLINK.txt
sed -i 's/N10/0/g' $DIR/CDHCall_ForPLINK.txt
#This drops the original column names
sed -i '1d' $DIR/CDHCall_ForPLINK.txt

#This adds a full column of 0 as the first column for Family ID
awk 'BEGIN{OFS=FS=" "} {$1=$1 OFS $1}1' $DIR/CDHCall_ForPLINK.txt > $DIR/CDHCall_ForPLINK_new.txt
mv $DIR/CDHCall_ForPLINK_new.txt $DIR/CDHCall_ForPLINK.txt

#This subsets the files for each CDH
awk -F' ' '{print $1, $2, $3}' $DIR/CDHCall_ForPLINK.txt > $DIR/Ab10Call_ForPLINK.txt
awk -F' ' '{print $1, $2, $4}' $DIR/CDHCall_ForPLINK.txt > $DIR/K10L2Call_ForPLINK.txt
awk -F' ' '{print $1, $2, $5}' $DIR/CDHCall_ForPLINK.txt > $DIR/BChromCall_ForPLINK.txt
awk -F' ' '{print $1, $2, $6}' $DIR/CDHCall_ForPLINK.txt > $DIR/CopyBChromCall_ForPLINK.txt

#This adds the column names for PLINK
sed -i '1i FID IID PHENOTYPE' $DIR/Ab10Call_ForPLINK.txt
sed -i '1i FID IID PHENOTYPE' $DIR/K10L2Call_ForPLINK.txt
sed -i '1i FID IID PHENOTYPE' $DIR/BChromCall_ForPLINK.txt
sed -i '1i FID IID PHENOTYPE' $DIR/CopyBChromCall_ForPLINK.txt

#This reads in my VCF and outputs a ped file. Double ID tells it that there are underscores in my names. I am intentionally not includeing a phenotype yet 
plink --vcf $DIR/$NAME.vcf.gz --recode --double-id --make-bed --out $DIR/$NAME

#This filters samples with >10% missing data
plink --mind 0.1 --double-id --make-bed --file $DIR/$NAME --out $DIR/$NAME

#This calculates a PCA for the samples to account for population structure
plink --bfile $DIR/$NAME --pca 10 'header' 'tabs' 'var-wts' --out $DIR/Ab10_PCA

#This performs a logistic GWAS for Ab10. --1 indicates my phenotype are coded as 0 and 1
plink --bfile $DIR/$NAME --logistic --beta -covar $DIR/Ab10_PCA.eigenvec --pheno $DIR/Ab10Call_ForPLINK.txt --1 --allow-no-sex --out Ab10_results_PCA
