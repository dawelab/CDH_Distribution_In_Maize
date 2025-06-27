#!/usr/bin/bash
#SBATCH --partition=highmem_p 
#SBATCH -J MAFFT_Ab10Hap
#SBATCH --output MAFFT_Ab10Hap.out
#SBATCH --mem=600GB
#SBATCH --time=72:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=30
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load MAFFT/7.520-GCC-12.3.0-with-extensions

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"

#Align the sequences

mafft --thread 30 $DIR/Consensus/HigginsN10_Ab10Hap.fa > $DIR/HigginsN10_Ab10Hap.MAFFTaln.fasta
