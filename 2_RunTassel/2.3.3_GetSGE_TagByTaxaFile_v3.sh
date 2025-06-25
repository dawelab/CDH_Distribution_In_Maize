#This defines variables
DIR=""

#I made this file in excel this converts it to a unix file, without this there are \r 
dos2unix $DIR/Get_K10L2_Table.txt

#This pulls info from the array job
IT=$SLURM_ARRAY_TASK_ID
DIRNAME=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $1}' $DIR/Get_K10L2_Table.txt)
TAXA=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $2}' $DIR/Get_K10L2_Table.txt)
SGE=$(awk NR==${SLURM_ARRAY_TASK_ID}'{print $3}' $DIR/Get_K10L2_Table.txt)


echo "This is iteration Number:"
echo $IT
echo "It is working on file:"
echo $DIRNAME
echo "It is extracting this SGE:"
echo $SGE

#This checks the length of the number of entries in the bed file
echo "##### Original Bed File Has This Number of Entries"
wc -l $DIR/$DIRNAME/$SGE"_TagID.bed"

#This isolates only the tag sequence from each line
cut -f 4 $DIR/$DIRNAME/$SGE"_TagID.bed" | cut -d "=" -f 2 > $DIR/$DIRNAME/$SGE"_temp1.bed"

#This checks the length of the number of entries in the bed file
echo "##### List Of Bed File Tags Has This Number Of Entries"
wc -l $DIR/$DIRNAME/$SGE"_temp1.bed"

#This checks the length of the original Taxa file
echo "##### Original Taxa File Has This Number of Entries"
wc -l $DIR/$DIRNAME/$TAXA".txt"

#This checks the length of the tag list from the Taxa file
echo "##### Taxa Tag List Has This Number of Entries"
wc -l $DIR/$DIRNAME/$TAXA"_TagsOnly.txt"

#This locates the index in the list of tags from the BIG tag by taxa file. The list of indexes comes from the BIG tag by taxa file so the index can still be used. 
#grep -n returns the index, grep -w returns only "whole word matches".
while read line
do
grep -w -n $line $DIR/$DIRNAME/$TAXA"_TagsOnly.txt" | cut -d ":" -f 1 >> $DIR/$DIRNAME/$SGE"_temp2.bed"
done < $DIR/$DIRNAME/$SGE"_temp1.bed"

#This checks the length of the SGE TAXA index list
echo "##### SGE Taxa File Index List Has This Number of Entries"
wc -l $DIR/$DIRNAME/$SGE"_temp2.bed"

#This isolates the header from the BIG tag by taxa file to add to each smaller file
head -n 1 $DIR/$DIRNAME/$TAXA".txt" > $DIR/$DIRNAME/$SGE"_temp3.bed"

#This uses the list of indexes and pulls the correct lines from the BIG tag by taxa file
awk 'FNR==NR{h[$1]; c++; next} FNR in h{print; if (!--c) exit}' $DIR/$DIRNAME/$SGE"_temp2.bed" $DIR/$DIRNAME/$TAXA".txt" >> $DIR/$DIRNAME/$SGE"_temp4.bed"

#This checks the length of the resulting SGE TAXA file
echo "##### SGE Taxa File Index List Has This Number of Entries"
wc -l $DIR/$DIRNAME/$SGE"_temp4.bed"

#This concatenates the header and the data
cat $DIR/$DIRNAME/$SGE"_temp3.bed" $DIR/$DIRNAME/$SGE"_temp4.bed" > $DIR/$DIRNAME/$TAXA.$SGE.txt

#This removes temporary files
rm $DIR/$DIRNAME/$SGE"_temp1.bed"
rm $DIR/$DIRNAME/$SGE"_temp2.bed"
rm $DIR/$DIRNAME/$SGE"_temp3.bed"
rm $DIR/$DIRNAME/$SGE"_temp4.bed"
