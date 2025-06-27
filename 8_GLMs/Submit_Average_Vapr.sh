#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Submit_Average_Vapr.sh
#SBATCH --output Submit_Average_Vapr.out
#SBATCH --mem=100GB
#SBATCH --time=48:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#Load the modules needed for sapelo2
module load R/4.4.1-foss-2022b

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Average_Vapr.R