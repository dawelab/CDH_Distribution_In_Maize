#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Submit_Ab10_KMeansModel_MixedControls_v2
#SBATCH --output Submit_Ab10_KMeansModel_MixedControls_v2.%A-%a.out
#SBATCH --mem=30GB
#SBATCH --time=10:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --array=1-25

#This defines the variables
IT=$SLURM_ARRAY_TASK_ID

#This writes out what file is working 
echo "working on file"
echo $IT

#Load the modules needed for sapelo2
module load R/4.3.1-foss-2022a

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/scripts/Ab10_KMeansModel_MixedControls_v2.R $IT
