#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J ID_SGE_Tags_v2
#SBATCH --output ID_SGE_Tags_v2.out
#SBATCH --mem=5GB
#SBATCH --time=5:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#This section separates the bed files into individual file

cd /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/
FILE="B73Ab10HiFiBChrom_Chunks.txt"
awk NR==1 $FILE > Ab10.bed
awk NR==2 $FILE > BChrom.bed

#This identifies only the tag field that overlaps with Ab10 or the B chromosome

module load BEDTools/2.30.0-GCC-10.2.0

cd /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2

bedtools intersect -wa -a BWAaln_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.s.bed -b Ab10.bed > Ab10_TagID.bed
bedtools intersect -wa -a BWAaln_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.s.bed -b BChrom.bed > BChrom_TagID.bed

awk '{print $4}' Ab10_TagID.bed > Ab10_TagIDOnly.bed
awk '{print $4}' BChrom_TagID.bed > BChrom_TagIDOnly.bed

#This gets the Tag only field from the giant tag by taxa file 
awk '{print $1}' Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2.txt > Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2_TagsOnly.txt

