#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J SelectLongestIsoform_K10L2
#SBATCH --output SelectLongestIsoform_K10L2.out
#SBATCH --mem=50GB
#SBATCH --time=2:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load SeqKit/0.16.1

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder"
PROT="/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/CI66_K10L2.v2.protein.fasta"

gunzip $DIR/$PROT.gz

while read i
do
seqkit faidx -r $PROT $i > $DIR/${i}_K10L2_protein.tmp.fasta
done < $DIR/CI66_K10L2.v2.protein.K10L2HapGeneNamesLongest.txt

cat $DIR/*_K10L2_protein.tmp.fasta > $DIR/CI66_K10L2.v2.protein.longest.fasta
rm $DIR/*_K10L2_protein.tmp.fasta