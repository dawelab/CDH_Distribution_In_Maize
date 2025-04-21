#!/bin/bash
#SBATCH --partition=batch
#SBATCH -J SelectK10L2Proteins
#SBATCH --output SelectK10L2Proteins.out
#SBATCH --mem=50G
#SBATCH --time=5:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load SeqKit/0.16.1 
module load AGAT/1.1.0

#Extract the protein sequence for K10L2 genes
#agat_sp_extract_sequences.pl -g /scratch/mjb51923/TRKIN_CRISPR/out_paper/Liftoff/CI66_K10L2_v1.gene.v2.sorted.gff3 -f /scratch/mjb51923/TRKIN_CRISPR/out_paper/RepeatMasker/CI66_K10L2_v1.fasta.masked -aa -o $DIR/CI66_K10L2.v2.protein.fasta

#Define the variables
DIR="/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder"
PROT=/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/CI66_K10L2.v2.protein.fasta

while read i
do
seqkit faidx -r $PROT $i > $DIR/${i}_K10L2Hap_protein.tmp.fasta
done < $DIR/CI66_K10L2_v1.gene.v2.K10L2hapGeneNames.txt

cat $DIR/*_K10L2Hap_protein.tmp.fasta > $DIR/CI66_K10L2.v2.protein.K10L2Hap.fasta
rm $DIR/*_K10L2Hap_protein.tmp.fasta