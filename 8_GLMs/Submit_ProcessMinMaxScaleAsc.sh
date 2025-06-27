#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Submit_ProcessMinMaxScaleAsc.sh
#SBATCH --output Submit_ProcessMinMaxScaleAsc.%A-%a.out
#SBATCH --mem=50GB
#SBATCH --time=1:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --array=1-24

#This defines the variables
IT=$SLURM_ARRAY_TASK_ID

#This generates a list of all the directories
#The "" line is there because it isn't reading the first line in for some reason
DIR=("" \
"Ab10_MaxEntStand4" \
"Ab10_MaxEntStand4_SSP126" \
"Ab10_MaxEntStand4_SSP245" \
"Ab10_MaxEntStand4_SSP370" \
"Ab10_MaxEntStand4_SSP585" \
"Ab10_MaxEntStand4_SSPCurrent" \
"BChrom_MaxEntStand4" \
"BChrom_MaxEntStand4_SSP126" \
"BChrom_MaxEntStand4_SSP245" \
"BChrom_MaxEntStand4_SSP370" \
"BChrom_MaxEntStand4_SSP585" \
"BChrom_MaxEntStand4_SSPCurrent" \
"K10L2_MaxEntStand4" \
"K10L2_MaxEntStand4_SSP126" \
"K10L2_MaxEntStand4_SSP245" \
"K10L2_MaxEntStand4_SSP370" \
"K10L2_MaxEntStand4_SSP585" \
"K10L2_MaxEntStand4_SSPCurrent" \
"Maize_MaxEntStand4" \
"Maize_MaxEntStand4_SSP126" \
"Maize_MaxEntStand4_SSP245" \
"Maize_MaxEntStand4_SSP370" \
"Maize_MaxEntStand4_SSP585" \
"Maize_MaxEntStand4_SSPCurrent")

TYPE=("" \
"Ab10" \
"Ab10" \
"Ab10" \
"Ab10" \
"Ab10" \
"Ab10" \
"BChrom" \
"BChrom" \
"BChrom" \
"BChrom" \
"BChrom" \
"BChrom" \
"K10L2" \
"K10L2" \
"K10L2" \
"K10L2" \
"K10L2" \
"K10L2" \
"Maize" \
"Maize" \
"Maize" \
"Maize" \
"Maize" \
"Maize")

MODEL=("" \
"Current_Full" \
"SSP126" \
"SSP245" \
"SSP370" \
"SSP585" \
"Current_SSP" \
"Current_Full" \
"SSP126" \
"SSP245" \
"SSP370" \
"SSP585" \
"Current_SSP" \
"Current_Full" \
"SSP126" \
"SSP245" \
"SSP370" \
"SSP585" \
"Current_SSP" \
"Current_Full" \
"SSP126" \
"SSP245" \
"SSP370" \
"SSP585" \
"Current_SSP")

#This reports the file that is being worked on
echo "working on file"
echo ${DIR[$IT]} 
echo ${TYPE[$IT]} 
echo ${MODEL[$IT]}

#Load the modules needed for sapelo2
module load R/4.4.1-foss-2022b

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/Ab10-Global-Survey/7_Map/ProcessMinMaxScaleAsc.R ${DIR[$IT]} ${TYPE[$IT]} ${MODEL[$IT]}