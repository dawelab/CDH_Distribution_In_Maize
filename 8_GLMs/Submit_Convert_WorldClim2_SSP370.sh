#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Submit_Convert_WorldClim2_SSP370
#SBATCH --output Submit_Convert_WorldClim2_SSP370.%A-%a.out
#SBATCH --mem=30GB
#SBATCH --time=1:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --array=1-55

#This defines the variables
IT=$SLURM_ARRAY_TASK_ID

#This pulls the file
FILE=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $1}' /scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/SSP370_TiffList.txt)

#This writes out what file is working 
echo "working on file"
echo $IT
echo $FILE

#Load the modules needed for sapelo2
module load  R/4.4.1-foss-2022b

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Convert_WorldClim2_SSP370.R $FILE