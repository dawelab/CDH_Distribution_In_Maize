#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Subset_VCF_Chr_Mo17
#SBATCH --output Subset_VCF_Chr_Mo17.out
#SBATCH --mem=100GB
#SBATCH --time=2:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=24
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#Load the module
module load BCFtools/1.15.1-GCC-11.3.0

#Define the Variables
DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"
VCF="SGEOnly_Mo17"

#This zips the vcf file
#bgzip $DIR/$VCF.vcf

#Index the new vcf file
#bcftools index $DIR/$VCF".vcf.gz"

#Subset the SNPs to only include those in the 1-9 range
#bcftools view  --threads 24 -r 1,2,3,4,5,6,7,8,9,10 -o $DIR/$VCF".chr.vcf.gz" $DIR/$VCF.vcf.gz

#Index the new vcf file
#bcftools index $DIR/$VCF".chr.vcf.gz"

#Subset to select only alleles with at least a 5% minor allele frequency, a depth of 2, and a genotype quality of 60, and have sample missingness of less than 25%
#bcftools filter --threads 24 -i 'FORMAT/DP > 3  && FORMAT/DP < 20 && MAF[0] > 0.05 && GQ > 60 && F_MISSING<0.25' -o $DIR/${VCF}.chr.filt.vcf.gz $DIR/${VCF}.chr.vcf.gz

#Index the new vcf file
bcftools index $DIR/${VCF}.chr.filt.vcf.gz

#Subset the SNPs to only the landraces
bcftools view  --threads 24 -o $DIR/$VCF".chr.filt.landrace.vcf.gz" $DIR/${VCF}.chr.filt.vcf.gz
