#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Explore_SNP_Number
#SBATCH --output Explore_SNP_Number.out
#SBATCH --mem=100GB
#SBATCH --time=24:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=24
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END


module load BCFtools/1.15.1-GCC-11.3.0

bcftools stats --threads 24 /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/GWAS_Mo17_Ab10.chr.filt2.vcf.gz