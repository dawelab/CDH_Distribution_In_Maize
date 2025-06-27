module load SeqKit/0.16.1
module load SAMtools/1.18-GCC-12.3.0

#Append the B chromosome to the complete Mo17 genome
seqkit faidx Zm-B73_B_CHROMOSOME-MBSC-1.0.fa chrB > ChrB_temp.fa
cat /scratch/mjb51923/ref_genomes/Zm-Mo17-REFERENCE-CAU-2.0.fa ChrB_temp.fa > Zm-Mo17-REFERENCE-CAU-2.0_chrB.fa

#Isolate the Ab10 and K10L2 haplotype
seqkit faidx -r /scratch/mjb51923/ref_genomes/Ab10_HiFi_v2_corrected.fa chr10 > Ab10_chr10.fa
#This bed file contains the single line chr10:141115174-195055488
samtools faidx -r Ab10_HiFi_v2_corrected_Ab10Hap.bed Ab10_chr10.fa > Ab10chr10_Ab10Hap.fa

#This bed file contains the single line K10L2:2730186-31891546
seqkit faidx -r /scratch/mjb51923/TRKIN_CRISPR/out_paper/RepeatMasker/CI66_K10L2_v1.fasta K10L2 > CI66_K10L2chr.fa
samtools faidx -r CI66_K10L2.bed CI66_K10L2chr.fa > CI66_K10L2chr_K10L2Hap.fa

#Truncate the Mo17 chromosome 10 at the colored 1 gene to avoid having duplicate shared regions
#The colored 1 gene in Mo17 is at 138838780 
seqkit faidx -r /scratch/mjb51923/ref_genomes/Zm-Mo17-REFERENCE-CAU-2.0.fa chr1 > Mo17_temp1.fa
#This pulls chr1 and 10 this line selects only chr1
head -n 5122265 Mo17_temp1.fa > Mo17_temp1_clip.fa
seqkit faidx -r /scratch/mjb51923/ref_genomes/Zm-Mo17-REFERENCE-CAU-2.0.fa chr2 > Mo17_temp2.fa
seqkit faidx -r /scratch/mjb51923/ref_genomes/Zm-Mo17-REFERENCE-CAU-2.0.fa chr3 > Mo17_temp3.fa
seqkit faidx -r /scratch/mjb51923/ref_genomes/Zm-Mo17-REFERENCE-CAU-2.0.fa chr4 > Mo17_temp4.fa
seqkit faidx -r /scratch/mjb51923/ref_genomes/Zm-Mo17-REFERENCE-CAU-2.0.fa chr5 > Mo17_temp5.fa
seqkit faidx -r /scratch/mjb51923/ref_genomes/Zm-Mo17-REFERENCE-CAU-2.0.fa chr6 > Mo17_temp6.fa
seqkit faidx -r /scratch/mjb51923/ref_genomes/Zm-Mo17-REFERENCE-CAU-2.0.fa chr7 > Mo17_temp7.fa
seqkit faidx -r /scratch/mjb51923/ref_genomes/Zm-Mo17-REFERENCE-CAU-2.0.fa chr8 > Mo17_temp8.fa
seqkit faidx -r /scratch/mjb51923/ref_genomes/Zm-Mo17-REFERENCE-CAU-2.0.fa chr9 > Mo17_temp9.fa
seqkit faidx -r /scratch/mjb51923/ref_genomes/Zm-Mo17-REFERENCE-CAU-2.0.fa chr10 > Mo17_temp10.fa
#The bed file contains the single line chr10:1-138838780
samtools faidx -r Mo17_NoChr10Shared.bed Mo17_temp10.fa > Mo17_temp10_cut.fa

#Combine all chromosomes
cat Mo17_temp1_clip.fa Mo17_temp2.fa Mo17_temp3.fa Mo17_temp4.fa Mo17_temp5.fa Mo17_temp6.fa Mo17_temp7.fa Mo17_temp8.fa Mo17_temp9.fa Mo17_temp10_cut.fa > Zm-Mo17-REFERENCE-CAU-2.0_NoChr10Shared.fa

#This adds Ab10 to the edited Mo17 genome so there is only one shared region of chr10
cat Zm-Mo17-REFERENCE-CAU-2.0_NoChr10Shared.fa Ab10chr10_Ab10Hap.fa > Zm-Mo17-REFERENCE-CAU-2.0_Ab10.fa
cat Zm-Mo17-REFERENCE-CAU-2.0_NoChr10Shared.fa CI66_K10L2chr_K10L2Hap.fa > Zm-Mo17-REFERENCE-CAU-2.0_K10L2.fa
