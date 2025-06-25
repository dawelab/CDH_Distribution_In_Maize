#This defines variables
DIR=""

#This pulls info from the array job
IT=$SLURM_ARRAY_TASK_ID
#This file is available under 2.2.2
DIRNAME=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $1}' $DIRGet_SGE_Table.txt)
TAXA=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $2}' $DIR/Get_SGE_Table.txt)

#This writes out what file is working 
echo "This is iteration Number"
echo $IT
echo "working on file"
echo $DIRNAME

join $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.15.txt $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.14.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.1
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.1 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.13.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.2
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.2 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.12.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.3
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.3 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.11.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.4
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.4 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.10.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.5
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.5 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.9.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.6
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.6 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.8.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.7
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.7 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.7.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.8
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.8 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.6.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.9
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.9 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.5.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.10
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.10 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.4.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.11
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.11 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.3.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.12
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.12 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.2.txtcd  > $DIR/$DIRNAME/Subset_TaxaFiles/temp.13
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.13 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Sum.1.txt > $DIR/$DIRNAME/$TAXA.Sub1Perc.Sum.txt

rm $DIR/$DIRNAME/Subset_TaxaFiles/temp.*
