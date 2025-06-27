#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J SplitVCF_ForBeagle_BChrom
#SBATCH --output SplitVCF_ForBeagle_BChrom.out
#SBATCH --mem=100GB
#SBATCH --time=24:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=24
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load BCFtools/1.15.1-GCC-11.3.0

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"
VCF=GWAS_Mo17_Bchr.chr.landrace.filt

for chr in {1..10}; do
#This subsets to just one chromosome
  bcftools view --threads 24 -r $chr $DIR/$VCF.vcf.gz -Oz -o $DIR/${VCF}_chr${chr}.vcf.gz
  tabix -p vcf $DIR/${VCF}_chr${chr}.vcf.gz
done
