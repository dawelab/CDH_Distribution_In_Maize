#This defines variables
IT=$SLURM_ARRAY_TASK_ID

echo "File Number"
echo $IT


#Load the modules needed for sapelo2
module load R/4.3.1-foss-2022a

Rscript --vanilla SCG_DataFiltering_RNLineMerging_v2.R $IT
