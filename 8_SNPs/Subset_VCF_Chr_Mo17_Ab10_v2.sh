#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Subset_VCF_Chr_Mo17_Ab10_v2
#SBATCH --output Subset_VCF_Chr_Mo17_Ab10_v2.out
#SBATCH --mem=100GB
#SBATCH --time=4:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=24
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#Load the module
module load BCFtools/1.15.1-GCC-11.3.0

#Define the Variables
DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"
VCF="GWAS_Mo17_Ab10"
CDH="AB10"

#This zips the vcf file
#bgzip $DIR/$VCF.vcf

#Index the new vcf file
#bcftools index $DIR/$VCF".vcf.gz"

#Subset the SNPs to only include those in the 1-9 range
#bcftools view  --threads 24 -r 1,2,3,4,5,6,7,8,9,10 -o $DIR/$VCF".chr.vcf.gz" $DIR/$VCF.vcf.gz

#This subsets to just chromosome 1 for exploratory purposes
#bcftools view  --threads 24 -r 1 -o $DIR/$VCF".chr1.vcf.gz" $DIR/$VCF.vcf.gz

#Index the new vcf file
#bcftools index $DIR/$VCF".chr.vcf.gz"

#I need to define Quality score in the header
#This isolates the header
#bcftools view -h $DIR/${VCF}.chr.vcf.gz > $DIR/${CDH}_hdr.txt

#I manually added the Quality score to the header

#This adds the new header
#bcftools reheader -h $DIR/${CDH}_hdr.txt $DIR/${VCF}.chr.vcf.gz > $DIR/${VCF}.chrfix.vcf.gz

#Subset the SNPs to only the high coverage lines
#bcftools view --threads 24 -S $DIR/Landrace_Only.txt -o $DIR/$VCF".chr.landrace.vcf.gz" $DIR/${VCF}.chrfix.vcf.gz

#Subset to remove any sites with missing data in more than 60% of the samples, less than 3 reads, or a minor allele frequency below 0.05%
bcftools filter --threads 24 -i 'FORMAT/DP > 3 && MAF[0] > 0.005 && F_MISSING < 0.4' -o $DIR/${VCF}.chr.landrace.filt.v2.vcf.gz $DIR/$VCF".chr.landrace.vcf.gz"

#Index the new vcf file
bcftools index $DIR/${VCF}.chr.landrace.filt.v2.vcf.gz

#This selects only biallelic SNPs
bcftools view --threads 24 -m2 -M2 -v snps  -o $DIR/${VCF}.chr.landrace.filt2.v2.vcf.gz $DIR/${VCF}.chr.landrace.filt.v2.vcf.gz

#Index the new vcf file
tabix $DIR/${VCF}.chr.landrace.filt2.v2.vcf.gz

######################################################
#This section repeats the subsetting for the CDH as a control
#bcftools view  --threads 24 -r $CDH -o $DIR/$VCF'.'$CDH'.vcf.gz' $DIR/$VCF.vcf.gz

#This adds the new header
#bcftools reheader -h $DIR/${CDH}_hdr.txt $DIR/$VCF'.'$CDH'.vcf.gz' > $DIR/$VCF'.'$CDH'.fix.vcf.gz'

#Subset to remove any sites with missing data in more than 70% of the samples
#bcftools filter --threads 24 -i 'F_MISSING < 0.3' -o $DIR/$VCF"."$CDH".filt.vcf.gz" $DIR/$VCF'.'$CDH'.fix.vcf.gz'

#Index the new vcf file
#bcftools index $DIR/$VCF"."$CDH".filt.vcf.gz"

#This selects only biallelic SNPs
#bcftools view --threads 24 -m2 -M2 -v snps  -o $DIR/$VCF"."$CDH".filt2.vcf.gz" $DIR/$VCF"."$CDH".filt.vcf.gz"