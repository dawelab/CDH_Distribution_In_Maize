#load module 
module load BEDTools/2.30.0-GCC-10.2.0

#define variables
DIR=""

#This section separates the bed files into individual file
cd $DIR
FILE="B73Ab10HiFiBChrom_Chunks.txt"
awk NR==1 $FILE > Ab10.bed
awk NR==2 $FILE > BChrom.bed

#This identifies only the tag field that overlaps with Ab10 or the B chromosome
#This bed file is output by TASSEL
bedtools intersect -wa -a BWAaln_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.s.bed -b Ab10.bed > Ab10_TagID.bed
bedtools intersect -wa -a BWAaln_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.s.bed -b BChrom.bed > BChrom_TagID.bed

awk '{print $4}' Ab10_TagID.bed > Ab10_TagIDOnly.bed
awk '{print $4}' BChrom_TagID.bed > BChrom_TagIDOnly.bed

#This gets the Tag only field from the large tag by taxa file 
awk '{print $1}' Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.txt > Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2_TagsOnly.txt

