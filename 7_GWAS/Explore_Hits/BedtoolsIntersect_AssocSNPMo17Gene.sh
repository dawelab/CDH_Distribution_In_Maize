module load BEDTools/2.31.0-GCC-12.3.0

DIR=/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2

awk '{print "chr"$1, $3, $3, $2, $4, $5, $6, $7, $8, $9}'  /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/Ab10_SigAssocSNP.txt| tail -n +2 | awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' > $DIR/Ab10_SigAssocSNP.bed

bedtools intersect -loj -a $DIR/Ab10_SigAssocSNP.bed -b /scratch/mjb51923/annotations/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.gff3 > $DIR/Ab10SigAssocSNP_Mo17Genes.bed


#Zm00014ba104130
#Zm00014ba176840


awk '{print "chr"$1, $3, $3, $2, $4, $5, $6, $7, $8, $9}'  /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChr_SigAssocSNP.txt| tail -n +2 | awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' > $DIR/BChr_SigAssocSNP.bed

bedtools intersect -loj -a $DIR/BChr_SigAssocSNP.bed -b /scratch/mjb51923/annotations/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.gff3 > $DIR/BChrSigAssocSNP_Mo17Genes.bed

#Zm00014ba033120
#Zm00014ba155010
#Zm00014ba218650
#Zm00014ba218800
#Zm00014ba255100
#Zm00014ba387780

awk '{print "chr"$1, $3, $3, $2, $4, $5, $6, $7, $8, $9}'  /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/K10L2_SigAssocSNP.txt| tail -n +2 | awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' > $DIR/K10L2_SigAssocSNP.bed

bedtools intersect -loj -a $DIR/K10L2_SigAssocSNP.bed -b /scratch/mjb51923/annotations/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.gff3 > $DIR/K10L2SigAssocSNP_Mo17Genes.bed
#Zm00014ba384340

module load SeqKit/0.16.1

PROT="Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.protein.fa"

for i in Zm00014ba104130 \
Zm00014ba176840 \
Zm00014ba033120 \
Zm00014ba155010 \
Zm00014ba218650 \
Zm00014ba218800 \
Zm00014ba384340 \
Zm00014ba255100 \
Zm00014ba387780
do
seqkit faidx -r $PROT $i > $DIR/${i}_CDHAssoc_protein.tmp.fasta
done 

cat *_CDHAssoc_protein.tmp.fasta > CDHAssoc_protein.fasta
rm *_CDHAssoc_protein.tmp.fasta
