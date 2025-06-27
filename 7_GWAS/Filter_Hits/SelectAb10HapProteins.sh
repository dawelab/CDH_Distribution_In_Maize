#!/bin/bash
#SBATCH --partition=batch
#SBATCH -J SelectAb10HapProteins
#SBATCH --output SelectAb10HapProteins.out
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=END

module load SeqKit/0.16.1
#module load AGAT/1.1.0

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder"

#Extract the protein sequence for K10L2 genes
#agat_sp_extract_sequences.pl -g /scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/Ab10_HiFi_v2_corrected.gene.v2.gff3 -f /scratch/mjb51923/TRKIN_CRISPR/out_paper/RepeatMasker/Ab10_HiFi_v2_corrected.fa.masked -aa -o $DIR/Ab10_HiFi_v2_corrected.gene.v2.proteins.fasta

PROT="/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/Ab10_HiFi_v2_corrected.gene.v2.proteins.fasta"

while read i
do
seqkit faidx -r $PROT $i > $DIR/${i}_Ab10Hap_protein.tmp.fasta
done < $DIR/Ab10_HiFi_v2_corrected.gene.v2.Ab10hapGeneNames.txt

cat $DIR/*_Ab10Hap_protein.tmp.fasta > $DIR/Ab10_HiFi_v2_corrected.gene.v2.protein.Ab10Hap.fa
rm $DIR/*_Ab10Hap_protein.tmp.fasta