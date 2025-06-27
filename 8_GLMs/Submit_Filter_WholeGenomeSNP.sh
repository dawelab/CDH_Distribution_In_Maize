#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Submit_FilterWholeGenomeSNP.sh
#SBATCH --output Submit_FilterWholeGenomeSNP.out
#SBATCH --mem=150GB
#SBATCH --time=5:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#Load the modules needed for sapelo2
module load R/4.4.1-foss-2022b

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Filter_WholeGenomeSNP.R 