#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J SelectLongestIsoform_Mo17
#SBATCH --output SelectLongestIsoform_Mo17.out
#SBATCH --mem=50GB
#SBATCH --time=2:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END


module load SeqKit/0.16.1

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder"
PROT="/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.protein.fa"

gunzip $DIR/$PROT.gz

while read i
do
seqkit faidx -r $PROT $i > $DIR/${i}_Mo17_protein.tmp.fasta
done < $DIR/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.protein.GeneNamesLongest.txt

cat $DIR/*_Mo17_protein.tmp.fasta > $DIR/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.protein.Longest.fa
rm $DIR/*_Mo17_protein.tmp.fasta