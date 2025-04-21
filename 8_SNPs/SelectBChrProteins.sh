#!/bin/bash
#SBATCH --partition=batch
#SBATCH -J SelectBChrProteins
#SBATCH --output SelectBChrProteins.out
#SBATCH --mem=50G
#SBATCH --time=5:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load SeqKit/0.16.1 

#Define the variables
DIR="/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder"
PROT=/scratch/mjb51923/annotations/Zm-B73_B_CHROMOSOME-MBSC-1.0_Zm00044a.1.protein.fa

while read i
do
seqkit faidx -r $PROT $i > $DIR/${i}_BChr_protein.tmp.fasta
done < $DIR/Zm-B73_B_CHROMOSOME-MBSC-1.0_Zm00044a.1.BChrGeneNames.txt

cat $DIR/*_BChr_protein.tmp.fasta > $DIR/Zm-B73_B_CHROMOSOME-MBSC-1.0_Zm00044a.1.protein.BChr.fa
rm $DIR/*_BChr_protein.tmp.fasta