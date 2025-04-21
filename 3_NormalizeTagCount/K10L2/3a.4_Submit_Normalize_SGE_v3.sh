#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Submit_Normalize_SGE_v3
#SBATCH --output Submit_Normalize_SGE_v3.%A-%a.out
#SBATCH --mem=200GB
#SBATCH --time=100:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=12
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --array=1-2

#This defines variables
DIR="/scratch/mjb51923/Ab10_Global_Survey/out"

#This pulls info from the array job
IT=$SLURM_ARRAY_TASK_ID
DIRNAME=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $1}' $DIR/Get_K10L2_Table.txt)
NAME=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $2}' $DIR/Get_K10L2_Table.txt)
SGE=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $3}' $DIR/Get_K10L2_Table.txt)


#This writes out what file is working 
echo "This is iteration Number"
echo $IT
echo "working on file"
echo $DIRNAME
echo "normalizing this SGE"
echo $SGE

#Load the modules needed for sapelo2
module load R/4.3.1-foss-2022a

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/scripts/Normalize_SGE.R $DIR/$DIRNAME $NAME $SGE > /scratch/mjb51923/Ab10_Global_Survey/scripts/Normalize_SGE_$IT.Rout
