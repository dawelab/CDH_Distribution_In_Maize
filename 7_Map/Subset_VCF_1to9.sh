#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Subset_VCF_1to9
#SBATCH --output Subset_VCF_1to9.out
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
VCF="CoordinateOnly_Ab10BChrom"

#Subset the SNPs to only include those in the 1-9 range
#bcftools view  --threads 24 -r 1,2,3,4,5,6,7,8,9  -S CoordinateLines.txt -o $DIR/$VCF".chr1to9.vcf.gz" $DIR/AllData_v5.vcf.gz

#Subset to select only alleles with at least a 1% minor allele frequency, a depth of 2, and a genotype quality of 20

#bcftools filter --threads 24 -i 'FORMAT/DP > 3  && FORMAT/DP < 20 && MAF[0] > 0.05 && GQ > 60' -o $DIR/${VCF}.chr1to9.filt2.vcf.gz $DIR/${VCF}.chr1to9.vcf.gz

#Load In the file created in the R script and subset it to only chromosome 1 so that things run a bit faster
bgzip $DIR/CoordinateOnly_Ab10BChrom.chr1to9.filt4.vcf
bcftools index $DIR/CoordinateOnly_Ab10BChrom.chr1to9.filt4.vcf.gz
bcftools view  --threads 24 -r 1 -o $DIR/$VCF".filt4.chr1.vcf.gz" $DIR/CoordinateOnly_Ab10BChrom.chr1to9.filt4.vcf.gz

