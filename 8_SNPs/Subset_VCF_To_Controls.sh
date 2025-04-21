#!/usr/bin/bash
#SBATCH --partition=batch 
#SBATCH -J Subset_VCF_To_HigControls
#SBATCH --output Subset_VCF_To_HigControls.out
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
bgzip AllData_v5.vcf
#index the vcf file
bcftools index AllData_v5.vcf.gz
#I need to add quality score to the header, this isolates the header
bcftools view -h AllData_v5.vcf.gz  > header.txt 

#I manually edited the header and updated it
bcftools reheader -h header.txt -o AllData_v5.fix.vcf.gz AllData_v5.vcf.gz

#This indexes the new file 
bcftools index AllData_v5.fix.vcf.gz

#Extract only the Ab10 haplotype and the Higgins control lines from the vcf file
bcftools view -R Ab10_VcfCompatible.bed -S Controls_Hig.txt -o AllData_v5.Ab10Hap.HigControls.vcf AllData_v5.fix.vcf.gz
bcftools view -R  Ab10_Shared_VcfCompatible.bed -S Controls_SharedTree.txt -o AllData_v5.Ab10Shared.ControlsST.vcf AllData_v5.fix.vcf.gz


