module load SAMtools/1.16.1-GCC-11.3.0

cd /scratch/mjb51923/ref_genomes

###############################

#This is downloads the T to T Mo17 genome assembly from the 2023 paper
#wget https://data.cyverse.org/dav-anon/iplant/home/laijs/Zm-Mo17-REFERENCE-CAU-2.0/Zm-Mo17-REFERENCE-CAU-2.0.fa.gz

#This unzips it
#gunzip Zm-Mo17-REFERENCE-CAU-2.0.fa.gz 

#This part finds the colored 1 (Zm00001eb429330 141187279..141196584 in B75 v5) in Mo17
#wget https://data.cyverse.org/dav-anon/iplant/home/laijs/Zm-Mo17-REFERENCE-CAU-2.0/Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.gff3
#maizegdb lists the Mo17 homolog for the B73 colored 1 gene as Zm00014ba089260
#grep Zm00014ba089260 Zm-Mo17-REFERENCE-CAU-2.0_Zm00014ba.gff3
#It is located at chr10:138839008-138846265

#This selects only the Mo17 chromosomes
#samtools faidx Zm-Mo17-REFERENCE-CAU-2.0.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 > Zm-Mo17-REFERENCE-CAU-2.0.noscaf.fa

#contig ptg000038l is the K10L2 haplotype

#This adds the Ab10 haplotype to the unedited Mo17 chromosomes. Contig ptg000038l is the K10L2 haplotype
cat Zm-Mo17-REFERENCE-CAU-2.0.noscaf.fa /scratch/mjb51923/CI66_Assembly/out/ptg000038l.fasta > Zm-Mo17-REFERENCE-CAU-2.0.noscaf.K10L2.fa
