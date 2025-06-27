#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Prep_TrkinDataForSRA2
#SBATCH --output Prep_TrkinDataForSRA2.out
#SBATCH --mem=3GB
#SBATCH --time=10:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

cd /scratch/mjb51923/raw_reads/NCBISRA_Submission/Submission_try2

#gzip  B73-Ab10_1.tar

split -b 90000m B73-Ab10_2.tar.gz B73-Ab10_2.tar.gz