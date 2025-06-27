#This defines the variables
IT=$SLURM_ARRAY_TASK_ID

#This writes out what file is working 
echo "working on file"
echo $IT

#Load the modules needed for sapelo2
module load R/4.3.1-foss-2022a

Rscript --vanilla K10L2_KMeansModel_MixedControls_v2.R $IT
