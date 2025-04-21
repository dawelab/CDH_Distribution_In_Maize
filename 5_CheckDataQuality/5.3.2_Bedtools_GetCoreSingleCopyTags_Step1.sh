#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Bedtools_GetCoreSingleCopyTags_Step1
#SBATCH --output Bedtools_GetCoreSingleCopyTags_Step1.out
#SBATCH --mem=400GB
#SBATCH --time=168:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#This loads the modules needed
module load BEDTools/2.30.0-GCC-10.2.0

#This defines the variables
DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"
ANNOT_DIR="/scratch/mjb51923/annotations"
BED="BWAaln_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.s"
TAXA="Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2"
CHUNKS="SingleCopyCoreGene_Chunks.txt"

#This checks the number of single copy core genes 
echo "There are this many single copy core genes"
wc -l $DIR/single_core_B73v5_genename.txt

#This selects the HiFiAb10 gene annotation with the same name as the B73 single copy core gene. It only selects the first occurance, liftoff created some fully duplicate lines. 
while read line
do
grep $line $ANNOT_DIR/Ab10_HiFi_v2_corrected.liftoff.gff3 | grep -m1 "gene" >> $DIR/Ab10_HiFi_v2_corrected.liftoff.CoreSingleCopy.gff3
done < $DIR/single_core_B73v5_genename.txt

#This checks the number of HiFiAb10 single copy core genes
echo "After selecting single copy core genes in the HifiAb10 genome there were this many "
wc -l $DIR/Ab10_HiFi_v2_corrected.liftoff.CoreSingleCopy.gff3

#This pulls tags that overlap with the single copy core genes
bedtools intersect -wa -a $DIR/$BED.bed -b $DIR/Ab10_HiFi_v2_corrected.liftoff.CoreSingleCopy.gff3 > $DIR/$BED.SingleCopyCoreGeneTags.s.bed

#This drops non unique lines 
sort $DIR/$BED.SingleCopyCoreGeneTags.s.bed | uniq > $DIR/$BED.SingleCopyCoreGeneTags.nodups.s.bed

#This extracts only the tag field and drops the tagSeq= string
awk '{print $4 }' $DIR/$BED.SingleCopyCoreGeneTags.nodups.s.bed | sed 's/tagSeq=//g' > $DIR/$BED.SingleCopyCoreGeneTags.TagOnly.txt

mkdir "SingleCopyCoreGeneChunks"

for i in {1..20}
do
START=$(awk NR==${i}'{print $1}' $DIR/$CHUNKS)
END=$(awk NR==${i}'{print $2}' $DIR/$CHUNKS)
for ((x=$START; x<=$END; x++)); { echo $x; } > $DIR/SingleCopyCoreGeneChunks/TempIndex.$i
awk 'FNR==NR{h[$1]; c++; next} FNR in h{print; if (!--c) exit}' $DIR/SingleCopyCoreGeneChunks/TempIndex.$i $DIR/$BED.SingleCopyCoreGeneTags.TagOnly.txt > $DIR/SingleCopyCoreGeneChunks/$BED.SingleCopyCoreGeneTags.TagOnly.$i.txt
done

rm $DIR/SingleCopyCoreGeneChunks/TempIndex.*
