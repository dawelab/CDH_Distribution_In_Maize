#!/bin/bash
#SBATCH --job-name=Submit_Ab10_DataFiltering_RNLineMerging
#SBATCH --output Submit_Ab10_DataFiltering_RNLineMerging.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=12
#SBATCH --mem=100gb
#SBATCH --time=48:00:00

#Load modules
module load R/4.3.1-foss-2022a

#Submit the script
R CMD BATCH /scratch/mjb51923/Ab10_Global_Survey/scripts/Ab10_DataFiltering_RNLineMerging.R
