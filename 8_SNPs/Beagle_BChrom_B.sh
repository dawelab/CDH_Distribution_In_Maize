#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Beagle_BChrom_B
#SBATCH --output Beagle_BChrom_B.out
#SBATCH --mem=100GB
#SBATCH --time=2:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=24
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#Load modules
module load Beagle/5.4.22Jul22.46e-Java-11
module load BCFtools/1.15.1-GCC-11.3.0

#Define variables
DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"
VCF="GWAS_Mo17_Bchr.B.landrace.filt"
CDH="B"

cd $DIR

#This performs imputation for beagle. This paper report using default paramets for maize https://pmc.ncbi.nlm.nih.gov/articles/PMC7532167/
java -Xmx90g -jar /apps/eb/Beagle/5.4.22Jul22.46e-Java-11/beagle.jar \
gt=$DIR/$VCF".vcf.gz" \
out=$DIR/$VCF".imputed" \
ne=1000 \
nthreads=24

tabix -p vcf $DIR/$VCF".imputed.vcf.gz"

