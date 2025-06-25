#This defines variables
DIR="/"

#This converts the excel filt to unix
dos2unix $DIR/Sum1Perc_Table.txt

#This pulls info from the array job
IT=$SLURM_ARRAY_TASK_ID
DIRNAME=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $1}' $DIR/AlignGBS_HiFiAb10Corrected_v2/Sum1Perc_Table.txt)
NUM=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $2}' $DIR/AlignGBS_HiFiAb10Corrected_v2/Sum1Perc_Table.txt)
NAME=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $3}' $DIR/AlignGBS_HiFiAb10Corrected_v2/Sum1Perc_Table.txt)

#This writes out what file is working 
echo "This is iteration Number"
echo $IT
echo "working on file"
echo $DIRNAME
echo "subset number"
echo $NUM

cd $DIR/$DIRNAME/Subset_TaxaFiles

module load R/4.3.1-foss-2022a

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/scripts/Sum_1Perc.R $DIR/$DIRNAME/Subset_TaxaFiles $NAME.$NUM.txt
