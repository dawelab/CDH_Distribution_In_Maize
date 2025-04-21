#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Submit_Filter_MaxEntAscFiles_CDHAssociated.sh
#SBATCH --output Submit_Filter_MaxEntAscFiles_CDHAssociated.%A-%a.out
#SBATCH --mem=50GB
#SBATCH --time=1:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --array=1-11

#This defines the variables
IT=$SLURM_ARRAY_TASK_ID

#This generates a list of all the directories
#The "" line is there because it isn't reading the first line in for some reason
DIR=("" \
"BChrom_MaxEntStand7" \
"BChrom_MaxEntStand7_SSP126" \
"BChrom_MaxEntStand7_SSP245" \
"BChrom_MaxEntStand7_SSP370" \
"BChrom_MaxEntStand7_SSP585" \
"K10L2_MaxEntStand7" \
"K10L2_MaxEntStand7_SSPCurrent" \
"K10L2_MaxEntStand7_SSP126" \
"K10L2_MaxEntStand7_SSP245" \
"K10L2_MaxEntStand7_SSP370" \
"K10L2_MaxEntStand7_SSP585")

TYPE=("" \
"BChrom" \
"BChrom" \
"BChrom" \
"BChrom" \
"BChrom" \
"K10L2" \
"K10L2" \
"K10L2" \
"K10L2" \
"K10L2"
"K10L2")

MODEL=("" \
"Current_Full" \
"SSP126" \
"SSP245" \
"SSP370" \
"SSP585" \
"Current_Full" \
"Current_SSP"
"SSP126" \
"SSP245" \
"SSP370" \
"SSP585")

#This reports the file that is being worked on
echo "working on file"
echo ${DIR[$IT]} 
echo ${TYPE[$IT]} 
echo ${MODEL[$IT]}

#Load the modules needed for sapelo2
module load R/4.4.1-foss-2022b

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/Filter_MaxEntAscFiles_CDHAssociated.R ${DIR[$IT]} ${TYPE[$IT]} ${MODEL[$IT]}