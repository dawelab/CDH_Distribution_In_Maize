#This defines the variables
DIR=""
ANNOT_DIR=""
#These files are from 2
BED="BWAaln_AllData_v5_v_B73-Ab10_BChrom.s"
TAXA="Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10_BChrom"
#This file is available in this repo under 5.4.2
CHUNKS="SingleCopyCoreGene_Chunks.txt"

#This concatenates all the Single Copy Core Gene Subsets
LIST=$(ls $DIR/SingleCopyCoreGeneChunks/Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10_BChrom_temp2*)
cat $LIST > $DIR/SingleCopyCoreGeneChunks/Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10_BChrom_temp2.all.bed

#This isolates the header from the BIG tag by taxa file to add to each smaller file
head -n 1 $DIR/$TAXA".txt" > $DIR/$TAXA"_header.bed"

#This adds the header to the file 
cat $DIR/$TAXA"_header.bed" $DIR/SingleCopyCoreGeneChunks/Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10_BChrom_temp2.all.bed > $DIR/Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10_BChrom_SingleCopyCoreGenes.bed

#This checks the length of list of tags within single copy core genes 
echo "##### There are this many tags overlapping single copy core genes"
wc -l $DIR/$BED.SingleCopyCoreGeneTags.TagOnly.txt

#This checks the length of the taxa file subset to only tags within single copy core genes
echo "##### There are this many tags in the tag by taxa file subset down to single copy core genes"
wc -l $DIR/Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10_BChrom.SingleCopyCoreGenes.txt
