#!/usr/bin/bash
#SBATCH --partition=batch 
#SBATCH -J Kalign_3_Shared
#SBATCH --output Kalign_3_Shared.out
#SBATCH --mem=150GB
#SBATCH --time=72:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=30
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load Kalign/3.3.2-GCCcore-11.2.0

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"

#Align the sequences
#kalign -i $DIR/Consensus/HigginsN10_Shared3.fa -o $DIR/HigginsN10_Shared3.kaln.fasta
#kalign -i $DIR/Consensus/HigginsN10_Ab10Hap.fa -o $DIR/HigginsN10_Ab10Hap.kaln.fasta
kalign -i $DIR/Consensus/HigginsN10_Shared1.fa -o $DIR/HigginsN10_Shared1.kaln.fasta 
#kalign -i $DIR/Consensus/HigginsN10_Shared2.fa -o $DIR/HigginsN10_Shared2.kaln.fasta
