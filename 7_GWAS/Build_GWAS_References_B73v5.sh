module load SeqKit/0.16.1
module load SAMtools/1.18-GCC-12.3.0

cd /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2


seqkit faidx -r /scratch/mjb51923/ref_genomes/Zm-B73_B_CHROMOSOME-MBSC-1.0.fa chrB > ChrB_temp.fa

seqkit faidx -r /scratch/mjb51923/ref_genomes/Ab10_HiFi_v2_corrected.fa chr10 > Ab10_chr10.fa
samtools faidx -r Ab10_HiFi_v2_corrected_Ab10Hap.bed Ab10_chr10.fa > Ab10chr10_Ab10Hap.fa

seqkit faidx -r /scratch/mjb51923/TRKIN_CRISPR/out_paper/RepeatMasker/CI66_K10L2_v1.fasta K10L2 > CI66_K10L2chr.fa
samtools faidx -r CI66_K10L2.bed CI66_K10L2chr.fa > CI66_K10L2chr_K10L2Hap.fa

#The colored 1 gene in B73v3 is at 138838780 
#This pulls chr1 and 10 this line selects only chr1
seqkit faidx -r /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/B73_RefGen_v3.fa Chr1 dna:chromosome chromosome:AGPv3:1:1:301476924:1 > B73v3_temp1.fa
#This pulls chr1 and 10 this line selects only chr1
tail -n +2493873 B73v3_temp1.fa > B73v3_temp1_clip.fa
seqkit faidx -r /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/B73_RefGen_v3.fa Chr2 dna:chromosome chromosome:AGPv3:2:1:237917468:1 > B73v3_temp2.fa
seqkit faidx -r /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/B73_RefGen_v3.fa Chr3 dna:chromosome chromosome:AGPv3:3:1:232245527:1 > B73v3_temp3.fa
seqkit faidx -r /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/B73_RefGen_v3.fa Chr4 dna:chromosome chromosome:AGPv3:4:1:242062272:1 > B73v3_temp4.fa
seqkit faidx -r /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/B73_RefGen_v3.fa Chr5 dna:chromosome chromosome:AGPv3:5:1:217959525:1 > B73v3_temp5.fa
seqkit faidx -r /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/B73_RefGen_v3.fa Chr6 dna:chromosome chromosome:AGPv3:6:1:169407836:1 > B73v3_temp6.fa
seqkit faidx -r /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/B73_RefGen_v3.fa Chr7 dna:chromosome chromosome:AGPv3:7:1:176826311:1 > B73v3_temp7.fa
seqkit faidx -r /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/B73_RefGen_v3.fa Chr8 dna:chromosome chromosome:AGPv3:8:1:175377492:1 > B73v3_temp8.fa
seqkit faidx -r /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/B73_RefGen_v3.fa Chr9 dna:chromosome chromosome:AGPv3:9:1:157038028:1 > B73v3_temp9.fa
seqkit faidx -r /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/B73_RefGen_v3.fa Chr10 dna:chromosome chromosome:AGPv3:10:1:149632204:1 > B73v3_temp10.fa
samtools faidx -r B73v3_NoChr10Shared.bed B73v3_temp10.fa > B73v3_temp10_cut.fa

cat B73v3_temp1_clip.fa B73v3_temp2.fa B73v3_temp3.fa B73v3_temp4.fa B73v3_temp5.fa B73v3_temp6.fa B73v3_temp7.fa B73v3_temp8.fa B73v3_temp9.fa B73v3_temp10_cut.fa > B73_RefGen_v3_NoChr10Shared.fa
cat B73v3_temp1_clip.fa B73v3_temp2.fa B73v3_temp3.fa B73v3_temp4.fa B73v3_temp5.fa B73v3_temp6.fa B73v3_temp7.fa B73v3_temp8.fa B73v3_temp9.fa B73v3_temp10.fa > B73_RefGen_v3_Chr.fa

#This add the B to the chr only reference
cat /scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/B73_RefGen_v3_Chr.fa ChrB_temp.fa > B73_RefGen_v3_Chr_chrB.fa
#This adds Ab10 to the edited B73v3 genome so there is only one shared region of chr10
cat B73_RefGen_v3_NoChr10Shared.fa Ab10chr10_Ab10Hap.fa > B73_RefGen_v3_NoChr10Shared_Ab10.fa
cat B73_RefGen_v3_NoChr10Shared.fa CI66_K10L2chr_K10L2Hap.fa > B73_RefGen_v3_NoChr10Shared_K10L2.fa
