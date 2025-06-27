#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Plink_GWAS_Landrace
#SBATCH --output Plink_GWAS_Landrace.out
#SBATCH --mem=100GB
#SBATCH --time=2:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load PLINK/1.9b_6.21-x86_64

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"
NAME="SGEOnly_Mo17.chr.filt.landrace"
PHENO="/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Controls_Swarts_RomeroNavarro_Romay_Groups_Env_Ab10K10L2BChromCopyNumCompleteWholeGenomePC.csv"
BChromList="/scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/8_SNPs/BChrom_Positive_Only.txt"

cd $DIR

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

#This investigates missing data
plink --file $DIR/$NAME --missing

#This filters SNPs with more than 25% missing data and samples with more than 10% missing data
plink --file $DIR/$NAME --geno 0.25 --mind 0.90 --make-bed --out $DIR/$NAME

#This investigates MAF
plink --bfile $DIR/$NAME --freq --out MAF_check

#This filters SNPs with MAF < 0.05 This was already done
#plink --bfile $DIR/$NAME --maf 0.05 --make-bed --noweb --out $DIR/$NAME

#This filters SNPs that are not in Hardy-Weinberg equilibrium. I am choosing not to do this because I don't know what to expect for suppressors and such
#plink --bfile $DIR/$NAME --hwe 1e-6 --make-bed --out $DIR/$NAME

#This removes SNPs that are in linkage
plink --bfile $DIR/$NAME --indep-pairwise 100000 50 0.5 --make-bed --out $DIR/$NAME
plink --bfile $DIR/$NAME --exclude $DIR/$NAME.prune.out --make-bed --out $DIR/$NAME

#This calculates a PCA for the samples to account for population structure
plink --bfile $DIR/$NAME --allow-extra-chr --pca 10 'header' 'tabs' 'var-wts'

#This performs a logistic GWAS for Ab10. --1 indicates my phenotype are coded as 0 and 
plink --bfile $DIR/$NAME --logistic --beta -covar $DIR/plink.eigenvec --pheno $DIR/Ab10Call_ForPLINK.txt --1 --allow-no-sex --out Ab10_results_PCA_Landrace

#This performs a logistic GWAS for K10L2. --1 indicates my phenotype are coded as 0 and 1
plink --bfile $DIR/$NAME --logistic --beta -covar $DIR/plink.eigenvec --pheno $DIR/K10L2Call_ForPLINK.txt --1 --allow-no-sex --out K10L2_results_PCA_Landrace

#This performs a logistic GWAS for K10L2. --1 indicates my phenotype are coded as 0 and 1
plink --bfile $DIR/$NAME --logistic --beta -covar $DIR/plink.eigenvec --pheno $DIR/BChromCall_ForPLINK.txt --1 --allow-no-sex --out BChrom_results_PCA_Landrace

#This performs a linear GWAS for BChrom pseudo copy number on only the B chromosome lines. 
plink --bfile $DIR/$NAME --linear --beta --keep $BChromList -covar $DIR/plink.eigenvec --pheno $DIR/CopyBChromCall_ForPLINK.txt --allow-no-sex --out CopyBChrom_results_PCA_Landrace