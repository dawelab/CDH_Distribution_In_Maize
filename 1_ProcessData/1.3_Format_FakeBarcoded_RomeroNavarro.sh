#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Format_FakeBarcoded_RomeroNavarro
#SBATCH --output /scratch/mjb51923/Ab10_Global_Survey/scripts/Format_FakeBarcoded_RomeroNavarro.out
#SBATCH --mem=100000
#SBATCH --time=100:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks-per-node=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=END

#Load the modules needed for sapelo2

DIR="/scratch/mjb51923/raw_reads/GBS/Romero_Navarro_Fake_Barcoded"

#Load modules


#This lists all of the folders in the directory 
cd $DIR
list=$(ls *.fastq) 

#This loops through all the files, alters the file name and then gzips the file ls 
for i in $list
do
mv $i $i.txt
gzip $i.txt
done
