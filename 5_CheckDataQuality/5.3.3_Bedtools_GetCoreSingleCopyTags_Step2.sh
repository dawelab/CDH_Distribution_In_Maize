#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Bedtools_GetCoreSingleCopyTags_Step2
#SBATCH --output Bedtools_GetCoreSingleCopyTags_Step2.%A-%a.out
#SBATCH --mem=100GB
#SBATCH --time=100:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --array=1-20

#This loads the modules needed
module load BEDTools/2.30.0-GCC-10.2.0

#This defines the variables
IT=$SLURM_ARRAY_TASK_ID
DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"
ANNOT_DIR="/scratch/mjb51923/annotations"
BED="BWAaln_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.s"
TAXA="Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2"
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
