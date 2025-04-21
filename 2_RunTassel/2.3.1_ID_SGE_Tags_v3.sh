#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J ID_SGE_Tags_v3
#SBATCH --output ID_SGE_Tags_v3.out
#SBATCH --mem=5GB
#SBATCH --time=5:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

#This Defines Variables
DIR="/scratch/mjb51923/Ab10_Global_Survey/out"
DIR_NAME="AlignGBS_Mo17K10L2"

cd $DIR/$DIR_NAME

#This identifies only the tag field that overlaps with Ab10 or the B chromosome

module load BEDTools/2.30.0-GCC-10.2.0

bedtools intersect -wa -a BWAaln_AllData_v6_v_Mo17K10L2.s.bed -b K10L2.bed > K10L2_TagID.bed

awk '{print $4}' K10L2_TagID.bed > K10L2_TagIDOnly.bed

#This gets the Tag only field from the giant tag by taxa file 
awk '{print $1}' Tassel_TagTaxaDist_AllData_v6_v_Mo17K10L2.txt > Tassel_TagTaxaDist_AllData_v6_v_Mo17K10L2_TagsOnly.txt
