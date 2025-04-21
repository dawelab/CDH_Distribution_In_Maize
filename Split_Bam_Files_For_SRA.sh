#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH -J Split_Bam_Files_For_SRA_1
#SBATCH --output Split_Bam_Files_For_SRA_1
#SBATCH --mem=30GB
#SBATCH --time=48:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=30
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=END

#Load modules
module load SAMtools/0.1.20-GCC-11.3.0

cd /scratch/mjb51923/raw_reads/NCBISRA_Submission/Submission_try2

samtools view -@ 30 -H m64060_210621_022359.subreads.bam > m64060_210621_022359.subreads.bam_header
samtools view -@ 30 m64060_210621_022359.subreads.bam | split -l 1000000 - m64060_210621_022359.subreads.split.bam.