#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J OrthoFinder_Mo17_K10L2
#SBATCH --output OrthoFinder_Mo17_K10L2.out
#SBATCH --mem=100GB
#SBATCH --time=24:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=30
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END


#load the modules necessary from the cluster
module load OrthoFinder/2.5.5-foss-2023a

#Define Variables 
DIR="/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder"

mkdir $DIR/Mo17_K10L2_proteomes
cd $DIR/Mo17_K10L2_proteomes
cp $DIR/CI66_K10L2.v2.protein.longest.fasta ./
cp $DIR/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.protein.Longest.fa ./

#This folder has the proteomes selected for the longest isoforms of Mo17 and BChr
orthofinder -t 30 -f $DIR/Mo17_K10L2_proteomes