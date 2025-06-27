#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Submit_Reformat_EnvRasters
#SBATCH --output Submit_Reformat_EnvRasters.out
#SBATCH --mem=50GB
#SBATCH --time=1:30:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#Load the modules needed for sapelo2
module load  R/4.4.1-foss-2022b

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Reformat_EnvRasters.R