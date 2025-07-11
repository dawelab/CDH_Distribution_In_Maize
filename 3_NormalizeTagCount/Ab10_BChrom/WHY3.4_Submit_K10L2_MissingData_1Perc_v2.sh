#!/bin/bash
#SBATCH --job-name=Submit_K10L2_MissingData_1Perc_v2
#SBATCH --output Submit_K10L2_MissingData_1Perc_v2.%A-%a.out
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mjb51923@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=24:00:00
#SBATCH --array=1-15

#This defines variables
DIR="/scratch/mjb51923/Ab10_Global_Survey/out"

#This converts dos (windows) line endings to unix line endings
dos2unix $DIR/AlignGBS_K10L2Contigs/K10L2_Sum1Perc_Table.txt

#This pulls info from the array job
IT=$SLURM_ARRAY_TASK_ID
DIRNAME=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $1}' $DIR/AlignGBS_K10L2Contigs/K10L2_Sum1Perc_Table.txt)
NUM=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $2}' $DIR/AlignGBS_K10L2Contigs/K10L2_Sum1Perc_Table.txt)
NAME=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $3}' $DIR/AlignGBS_K10L2Contigs/K10L2_Sum1Perc_Table.txt)

#This writes out what file is working 
echo "This is iteration Number"
echo $IT
echo "working on file"
echo $DIRNAME
echo "subset number"
echo $NUM

cd $DIR/$DIRNAME/Subset_TaxaFiles

module load R/4.3.1-foss-2022a

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/scripts/MissingData_1Perc.R $DIR/$DIRNAME/Subset_TaxaFiles $NAME.$NUM.txt
