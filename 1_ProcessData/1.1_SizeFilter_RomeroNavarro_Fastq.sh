DIR=""
OUT_DIR=""

#Load modules
module load seqtk/1.3-GCC-8.3.0

#List all fastq files in dir
cd $DIR
list=$(ls *.fastq) 

#Initiate a text file that I will add values too in the loop
printf "Line\tOrig_Count\tFilt_Count\n" > $OUT_DIR/RomeroNavarro_LengthFiltReadNumber.txt

#Loop through all the files and select only reads 64bp or longer
for i in $list
do
seqtk seq -L64 $DIR/$i > $OUT_DIR/$i
READS=$(echo $(cat $DIR/$i|wc -l)/4|bc)
FILT_READS=$(echo $(cat $OUT_DIR/$i|wc -l)/4|bc)
printf "${i}\t${READS}\t${FILT_READS}\n" >> $OUT_DIR/RomeroNavarro_LengthFiltReadNumber.txt
done
