#load the modules necessary from the cluster
module load BLAST+/2.13.0-gompi-2022a

#Define Variables, variables in all caps 
DIR="/scratch/mjb51923/annotations"
REF="/scratch/mjb51923/ref_genomes/Zm-B73_B_CHROMOSOME-MBSC-1.0.fa"

#make blast database for the B73Ab10 reference genome
makeblastdb -in $REF -parse_seqids -dbtype nucl

#blast trkin against the scaffolds 
blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query $DIR/"CL-repeat.fasta" -out $DIR/BLAST_CLrepeat_v_BChrom.out
blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query $DIR/"StarkB.fasta" -out $DIR/BLAST_StarkB_v_BChrom.out
