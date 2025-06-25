#This defines variables
DIR=""

#This pulls info from the array job
IT=$SLURM_ARRAY_TASK_ID
#The Get_SGE_Table.txt is available in 2.2.2
DIRNAME=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $1}' $DIR/Get_SGE_Table.txt)
TAXA=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $2}' $DIR/Get_SGE_Table.txt)

#This writes out what file is working 
echo "This is iteration Number"
echo $IT
echo "working in directory"
echo $DIRNAME
echo "on this file"
exho $TAXA

#This counts the number of rows in the Taxa file
NROW=$(wc -l $DIR/$DIRNAME/$TAXA.txt | awk '{print $1}')

#This calculates 1% of all the rows in the Taxa file
SAMP=$(echo "$NROW*.01" | bc)

#This extracts only the header of the Taxa file
head -n 1 $DIR/$DIRNAME/$TAXA.txt > $DIR/$DIRNAME/temp1.txt

#This makes a directory for all the tag index lists
mkdir $DIR/$DIRNAME/Subset_IndexLists
mkdir $DIR/$DIRNAME/Subset_TaxaFiles

#This loop extracts 1 percent of the Taxa file 15 times to get an accurate estimate of 1%
i=1
while (($i<=15 ))
do
 echo $i
 #This randomly selects indexes for 1% of all tags in the TAXA file and saves the lists for reproduction if needed
  shuf -i 1-$NROW -n $SAMP > $DIR/$DIRNAME/Subset_IndexLists/list_1Perc_$i.num
  awk 'FNR==NR{h[$1]; c++; next} FNR in h{print; if (!--c) exit}' $DIR/$DIRNAME/Subset_IndexLists/list_1Perc_$i.num $DIR/$DIRNAME/$TAXA.txt >> $DIR/$DIRNAME/Subset_TaxaFiles/temp2.$i.txt
  #This adds the header to the file
  cat $DIR/$DIRNAME/temp1.txt $DIR/$DIRNAME/Subset_TaxaFiles/temp2.$i.txt > $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.$i.txt
  #This adds one to the counter
  (( i=$i+1 ))
done

rm $DIR/$DIRNAME/temp1.txt
rm $DIR/$DIRNAME/Subset_TaxaFiles/temp2.*.txt

ls $DIR/$DIRNAME/Subset_TaxaFiles/
