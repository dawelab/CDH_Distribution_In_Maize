module load BEDTools/2.31.0-GCC-12.3.0

DIR=/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2

awk 'NR > 1 {print $6, $7, $8}' /scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/Mo17_Ab10_proteomes/Mo17GenesOrthologousToAb10.txt | tail -n +2 | awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' > $DIR/Mo17GenesOrthologousToAb10.bed
awk 'NR > 1 {print "chr"$1, $3, $3, $2, $4, $5, $6, $7, $8, $9}' /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/Ab10_results_PCA.assoc.logistic | tail -n +2 | awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' > $DIR/Ab10_results_PCA.assoc.logistic.bed

bedtools intersect -loj -a $DIR/Ab10_results_PCA.assoc.logistic.bed -b $DIR/Mo17GenesOrthologousToAb10.bed > $DIR/Ab10_results_PCA.assoc.logistic.Ab10Ortho.bed 



awk 'NR > 1 {print $6, $7, $8}' /scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/Mo17_K10L2_proteomes/Mo17GenesOrthologousToK10L2.txt | tail -n +2 | awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' > $DIR/Mo17GenesOrthologousToK10L2.bed
awk 'NR > 1 {print "chr"$1, $3, $3, $2, $4, $5, $6, $7, $8, $9}' /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/K10L2_results_PCA.assoc.logistic | tail -n +2 | awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' > $DIR/K10L2_results_PCA.assoc.logistic.bed

bedtools intersect -loj -a $DIR/K10L2_results_PCA.assoc.logistic.bed -b $DIR/Mo17GenesOrthologousToK10L2.bed > $DIR/K10L2_results_PCA.assoc.logistic.K10L2Ortho.bed 



awk 'NR > 1 {print $6, $7, $8}' /scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/Mo17_BChr_proteomes/Mo17GenesOrthologousToBChr.txt | tail -n +2 | awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' > $DIR/Mo17GenesOrthologousToBChr.bed
awk 'NR > 1 {print "chr"$1, $3, $3, $2, $4, $5, $6, $7, $8, $9}' /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/BChrom_results_PCA.assoc.logistic | tail -n +2 | awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' > $DIR/BChrom_results_PCA.assoc.logistic.bed

bedtools intersect -loj -a $DIR/BChrom_results_PCA.assoc.logistic.bed -b $DIR/Mo17GenesOrthologousToBChr.bed > $DIR/BChrom_results_PCA.assoc.logistic.BChromOrtho.bed 
