#!/usr/bin/bash
#SBATCH --partition=batch 
#SBATCH -J Subset_VCF_To_Ab10Pos
#SBATCH --output Subset_VCF_To_Ab10Pos.out
#SBATCH --mem=100GB
#SBATCH --time=5:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#Load the module
module load BCFtools/1.15.1-GCC-11.3.0

#Enter the working directory
cd /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2
#Zip the vcf file
bgzip SGEOnly_Ab10BChrom.vcf
#index the vcf file
bcftools index SGEOnly_Ab10BChrom.vcf.gz

#I need to add quality score to the header, this isolates the header
bcftools view -h SGEOnly_Ab10BChrom.vcf.gz > header.txt 

#I manually edited the header and updated it
bcftools reheader -h header.txt -o SGEOnly_Ab10BChrom.fix.vcf.gz SGEOnly_Ab10BChrom.vcf.gz

#This indexes the new file 
bcftools index SGEOnly_Ab10BChrom.fix.vcf.gz

#Extract only the Ab10 haplotype and the Higgins control lines from the vcf file
bcftools view -R Ab10_VcfCompatible.bed -S All_Ab10_Positive.txt -o SGEOnly_Ab10BChrom.Ab10Hap.AllAb10Pos.vcf.gz SGEOnly_Ab10BChrom.fix.vcf.gz


