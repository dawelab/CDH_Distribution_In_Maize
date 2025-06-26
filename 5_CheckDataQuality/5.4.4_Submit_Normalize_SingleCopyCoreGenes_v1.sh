#This defines variables
DIR=""
DIRNAME=""
NAME="Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10_BChrom"
IT=$SLURM_ARRAY_TASK_ID
SGE="SingleCopyCoreGenes"

#This writes out what file is working 
echo "working on file"
echo $DIRNAME
echo "normalizing this SGE"
echo $SGE

#Load the modules needed for sapelo2
module load R/4.3.1-foss-2022a

Rscript --vanilla Normalize_SingleCopyCoreGenes.R $DIR/$DIRNAME $NAME $IT $SGE
