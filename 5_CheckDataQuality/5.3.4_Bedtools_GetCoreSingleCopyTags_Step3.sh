#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Bedtools_GetCoreSingleCopyTags_Step3
#SBATCH --output Bedtools_GetCoreSingleCopyTags_Step3.out
#SBATCH --mem=5GB
#SBATCH --time=5:00:00
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

#This concatenates all the Single Copy Core Gene Subsets
LIST=$(ls $DIR/SingleCopyCoreGeneChunks/Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2_temp2*)
cat $LIST > $DIR/SingleCopyCoreGeneChunks/Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2_temp2.all.bed

#This isolates the header from the BIG tag by taxa file to add to each smaller file
head -n 1 $DIR/$TAXA".txt" > $DIR/$TAXA"_header.bed"

#This adds the header to the file 
cat $DIR/$TAXA"_header.bed" $DIR/SingleCopyCoreGeneChunks/Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2_temp2.all.bed > $DIR/Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2_SingleCopyCoreGenes.bed

#This checks the length of list of tags within single copy core genes 
echo "##### There are this many tags overlapping single copy core genes"
wc -l $DIR/$BED.SingleCopyCoreGeneTags.TagOnly.txt

#This checks the length of the taxa file subset to only tags within single copy core genes
echo "##### There are this many tags in the tag by taxa file subset down to single copy core genes"
wc -l $DIR/Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.SingleCopyCoreGenes.txt
