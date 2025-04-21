#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Submit_SCG_DataFiltering_RNLineMerging_v2
#SBATCH --output Submit_SCG_DataFiltering_RNLineMerging_v2.%A-%a.out
#SBATCH --mem=100GB
#SBATCH --time=24:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=12
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --array=1-20

#This defines variables
IT=$SLURM_ARRAY_TASK_ID

echo "File Number"
echo $IT


#Load the modules needed for sapelo2
module load R/4.3.1-foss-2022a

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/scripts/SCG_DataFiltering_RNLineMerging_v2.R $IT
