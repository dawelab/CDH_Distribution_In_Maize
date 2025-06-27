#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Plink_GWAS_BChr
#SBATCH --output Plink_GWAS_BChr.out
#SBATCH --mem=50GB
#SBATCH --time=2:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load PLINK/1.9b_6.21-x86_64

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"
PHENO="/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/8_SNPs/BChr_ContExperimental_CopyNumber.csv"

#This concatenates the imputed files
cd $DIR

bcftools concat -o $DIR/GWAS_Mo17_Bchr.chr.landrace.filt.imputed.filt.vcf.gz $DIR/GWAS_Mo17_Bchr.chr.landrace.filt_chr*.imputed.vcf.gz

NAME="GWAS_Mo17_Bchr.chr.landrace.filt.imputed.filt"

#This extracts the CDH call columns from the Pheno file
awk -F, '{print $1, $2}' $PHENO > $DIR/CopyNumCall_ForPLINK.txt

#This drops the original column names
sed -i '1d' $DIR/CopyNumCall_ForPLINK.txt

#This adds a full column of 0 as the first column for Family ID
awk 'BEGIN{OFS=FS=" "} {$1=$1 OFS $1}1' $DIR/CopyNumCall_ForPLINK.txt > $DIR/CopyNumCall_ForPLINK_new.txt
mv $DIR/CopyNumCall_ForPLINK_new.txt $DIR/CopyNumCall_ForPLINK.txt

#This adds the column names for PLINK
sed -i '1i FID IID PHENOTYPE' $DIR/CopyNumCall_ForPLINK.txt

#This reads in my VCF and outputs a ped file. Double ID tells it that there are underscores in my names. I am intentionally not includeing a phenotype yet 
plink --vcf $DIR/$NAME.vcf.gz --recode --double-id --make-bed --out $DIR/$NAME

#This filters samples with >10% missing data
plink --mind 0.1 --double-id --make-bed --file $DIR/$NAME --out $DIR/$NAME

#This performs a logistic GWAS for K10L2. --1 indicates my phenotype are coded as 0 and 1
plink --bfile $DIR/$NAME --linear --beta -covar $DIR/BChr_PCA.eigenvec --pheno $DIR/CopyNumCall_ForPLINK.txt --allow-no-sex --out BChromCopyNum_results_PCA
