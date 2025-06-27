#!/bin/bash
#!/usr/bin/bash
#SBATCH --partition=batch  
#SBATCH -J SubsetVCF1Perc.sh
#SBATCH --output SubsetVCF1Perc.out
#SBATCH --mem=100GB
#SBATCH --time=10:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load BCFtools/1.15.1-GCC-11.3.0

# Input and output VCF file names
INPUT_VCF=CoordinateOnly_Ab10BChrom.chr1to9.filt4.vcf.gz
OUTPUT_VCF=CoordinateOnly_Ab10BChrom.chr1to9.filt4.1Perc.vcf.gz

# Randomly subset 1% of variants using bcftools
bcftools view -Oz -r "$(bcftools view -H $INPUT_VCF | awk 'BEGIN {srand()} !/^#/ {if (rand() <= 0.01) print $1":"$2}' | tr '\n' ',')" $INPUT_VCF > $OUTPUT_VCF

# Index the output VCF file
bcftools index $OUTPUT_VCF

echo "Subset VCF file saved as: $OUTPUT_VCF"