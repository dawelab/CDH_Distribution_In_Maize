#This loads the modules needed
module load BEDTools/2.30.0-GCC-10.2.0

#This defines the variables
IT=$SLURM_ARRAY_TASK_ID
DIR=""
ANNOT_DIR=""
#These files are from 2
BED="BWAaln_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.s"
TAXA="Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2"
#This file is available in this repo under 5.4.3
CHUNKS="SingleCopyCoreGene_Chunks.txt"

#This checks the length of the SGE TAXA index subset list
echo "##### Taxa File Index Subset List Has This Number of Entries"
wc -l $DIR/SingleCopyCoreGeneChunks/$BED.SingleCopyCoreGeneTags.TagOnly.$IT.txt

#This locates the index in the list of tags from the BIG tag by taxa file. The list of indexes comes from the BIG tag by taxa file so the index can still be used. 
#grep -n returns the index, grep -w returns only "whole word matches".
while read line
do
grep -w -n $line $DIR/$TAXA"_TagsOnly.txt" | cut -d ":" -f 1 >> $DIR/SingleCopyCoreGeneChunks/$TAXA"_temp1_"$IT".bed"
done < $DIR/SingleCopyCoreGeneChunks/$BED.SingleCopyCoreGeneTags.TagOnly.$IT.txt

#This checks the length of the SGE TAXA index subset list
echo "##### Taxa File Subset Has This Number of Entries"
wc -l $DIR/SingleCopyCoreGeneChunks/$TAXA"_temp1_"$IT".bed"

#This uses the list of indexes and pulls the correct lines from the BIG tag by taxa file
awk -v FILE=$DIR/SingleCopyCoreGeneChunks/$TAXA"_temp1_"$IT".bed" 'BEGIN { while(getline < FILE) A[$1]=1 } A[NR]' < $DIR/$TAXA".txt" > $DIR/SingleCopyCoreGeneChunks/$TAXA"_temp2_"$IT".bed"

#This checks the length of the resulting SGE TAXA subset file
echo "##### Single Copy Core Gene Taxa File Has This Number of Entries"
wc -l $DIR/SingleCopyCoreGeneChunks/$TAXA"_temp2_"$IT".bed"
