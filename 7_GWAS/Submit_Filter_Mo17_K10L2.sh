#!/usr/bin/bash
#SBATCH --partition=highmem_p
#SBATCH -J Submit_Filter_Mo17_K10L2
#SBATCH --output Submit_Filter_Mo17_K10L2.out
#SBATCH --mem=400GB
#SBATCH --time=48:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#Load the modules needed for sapelo2
module load R/4.4.1-foss-2022b

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/8_SNPs/Filter_Mo17_K10L2.R