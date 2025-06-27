module load BLAST+/2.13.0-gompi-2022a
DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"
REF=/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/Zm-Mo17-REFERENCE-CAU-2.0_NoChr10Shared.fa

makeblastdb -in $REF -parse_seqids -dbtype nucl

blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query /scratch/mjb51923/annotations/KindrE9.fasta -out $DIR/BLAST_Kindr_v_Mo17.out

REF=/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/Zm-Mo17-REFERENCE-CAU-2.0_NoChr10Shared.fa

makeblastdb -in $REF -parse_seqids -dbtype nucl

blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query /scratch/mjb51923/annotations/KindrE9.fasta -out $DIR/BLAST_Kindr_v_Mo17.out
