
DIR="/path/to/barcode/faked/fastqfiles"

#This lists all of the files output by barcode faker
cd $DIR
list=$(ls *.fastq) 

#This loops through all the files, alters the file name and then gzips the file
for i in $list
do
mv $i $i.txt
gzip $i.txt
done
