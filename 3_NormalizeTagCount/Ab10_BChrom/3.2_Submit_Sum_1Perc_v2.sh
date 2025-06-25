#This defines variables
DIR="/"

#This converts the excel filt to unix
dos2unix $DIR/Sum1Perc_Table.txt

#This pulls info from the array job
IT=$SLURM_ARRAY_TASK_ID
#The Sum1Perc_Table file is provided in this directory under 3.2 as well
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

Rscript --vanilla Sum_1Perc.R $DIR/$DIRNAME/Subset_TaxaFiles $NAME.$NUM.txt
