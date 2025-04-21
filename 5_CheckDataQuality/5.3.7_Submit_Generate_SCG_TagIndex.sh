#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Submit_Generate_SCG_TagIndex
#SBATCH --output Submit_Generate_SCG_TagIndex.out
#SBATCH --mem=300GB
#SBATCH --time=24:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=12
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#Load the modules needed for sapelo2
module load R/4.3.1-foss-2022a

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/scripts/Generate_SCG_TagIndex.R
