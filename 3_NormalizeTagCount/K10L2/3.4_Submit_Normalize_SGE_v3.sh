#This pulls info from the array job
IT=$SLURM_ARRAY_TASK_ID
#This table is available under 2.3.2
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

Rscript --vanilla Normalize_SGE.R $DIR/$DIRNAME $NAME $SGE > Normalize_SGE_$IT.Rout
