#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J SizeFilter_RomeroNavarro_Fastq
#SBATCH --output /scratch/mjb51923/Ab10_Global_Survey/scripts/SizeFilter_RomeroNavarro_Fastq.out
#SBATCH --mem=100000
#SBATCH --time=150:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks-per-node=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=END

#Load the modules needed for sapelo2

DIR="/scratch/mjb51923/raw_reads/GBS/Romero_Navarro"
OUT_DIR="/scratch/mjb51923/raw_reads/GBS/Romero_Navarro_LengthFilt"

#Load modules
module load seqtk/1.3-GCC-8.3.0

#This lists all of the folders in the directory 
cd $DIR
list=$(ls *.fastq) 

#Initiate a text file that I will add values too in the loop
printf "Line\tOrig_Count\tFilt_Count\n" > $OUT_DIR/RomeroNavarro_LengthFiltReadNumber.txt

#This loops through all the files and size selects them  
for i in $list
do
seqtk seq -L64 $DIR/$i > $OUT_DIR/$i
READS=$(echo $(cat $DIR/$i|wc -l)/4|bc)
FILT_READS=$(echo $(cat $OUT_DIR/$i|wc -l)/4|bc)
printf "${i}\t${READS}\t${FILT_READS}\n" >> $OUT_DIR/RomeroNavarro_LengthFiltReadNumber.txt
done
