#This defines variables
DIR=""

#This pulls info from the array job
IT=$SLURM_ARRAY_TASK_ID
#This file is available in 3.2
DIRNAME=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $1}' $DIR/Sum1Perc_Table.txt)
NUM=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $2}' $DIR/Sum1Perc_Table.txt)
NAME=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $3}' $DIR/Sum1Perc_Table.txt)

#This writes out what file is working 
echo "This is iteration Number"
echo $IT
echo "working on file"
echo $DIRNAME
echo "subset number"
echo $NUM

cd $DIR/$DIRNAME/Subset_TaxaFiles

module load R/4.3.1-foss-2022a

Rscript --vanilla MissingData_1Perc.R $DIR/$DIRNAME/Subset_TaxaFiles $NAME.$NUM.txt
