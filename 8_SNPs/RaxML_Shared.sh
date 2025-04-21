#!/usr/bin/bash
#SBATCH --partition=batch 
#SBATCH -J RaxML_Shared
#SBATCH --output RaxML_Shared.out
#SBATCH --mem=50GB
#SBATCH --time=24:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=10
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load RAxML-NG/1.2.2-GCC-12.2.0

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"

#This command generates a ML tree using the GTR+G model in RAxML
#raxml-ng --check --all --msa $DIR/HigginsN10_Shared3.kaln.fasta --model "GTR+G" --seed 10 --prefix HigginsN10_Shared3 --threads 10
#raxml-ng --check --all --msa $DIR/HigginsN10_Shared2.kaln.fasta --model "GTR+G" --seed 10 --prefix HigginsN10_Shared2 --threads 10
raxml-ng --all --msa $DIR/HigginsN10_Shared1.kaln.fasta --model "GTR+G" --seed 10 --prefix HigginsN10_Shared1
#raxml-ng --check --all --msa $DIR/HigginsN10_Ab10Hap.kaln.fasta --model "GTR+G" --seed 10 --prefix HigginsN10_Ab10Hap --threads 10