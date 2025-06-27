#!/usr/bin/bash
#SBATCH --partition=highmem_p 
#SBATCH -J MUSCLE_Ab10Hap
#SBATCH --output MUSCLE_Ab10Hap.out
#SBATCH --mem=600GB
#SBATCH --time=72:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=30
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load MUSCLE/5.1.0-GCCcore-11.3.0

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"

#Align the sequences

muscle -super5 $DIR/Consensus/HigginsN10_Ab10Hap.fa -output $DIR/HigginsN10_Ab10Hap.kaln.fasta -threads 30
