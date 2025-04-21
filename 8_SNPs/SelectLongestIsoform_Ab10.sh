#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J SelectLongestIsoform_Ab10
#SBATCH --output SelectLongestIsoform_Ab10.out
#SBATCH --mem=50GB
#SBATCH --time=2:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load SeqKit/0.16.1

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder"
PROT="/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/Ab10_HiFi_v2_corrected.gene.v2.protein.Ab10Hap.fa"

while read i
do
i=${i//_NA/}
seqkit faidx -r $PROT $i > $DIR/${i}_Ab10_protein.tmp.fasta
done < $DIR/Ab10_HiFi_v2_corrected.gene.v2.Ab10hapGeneNamesLongest.txt

cat $DIR/*_Ab10_protein.tmp.fasta > $DIR/Ab10_HiFi_v2_corrected.gene.v2.protein.Ab10Hap.longest.fa
rm $DIR/*_Ab10_protein.tmp.fasta