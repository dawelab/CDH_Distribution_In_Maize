module load SAMtools/1.18-GCC-12.3.0

REF="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/Zm-Mo17-REFERENCE-CAU-2.0_chrB.fa"
DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"

#Write a line of code to split this 
	
Chr=chr3
START=$((5366738-64))
END=$((5366738+64))

samtools faidx $REF $Chr:$START-$END > $DIR/Mo17_${Chr}_${START}-${END}.fa


#load the modules necessary from the cluster
#module load BLAST+/2.13.0-gompi-2022a

#make blast database
#makeblastdb -in $REF -parse_seqids -dbtype nucl
#Blast various relevant features to the CI66_K10L2 genome and convert them to a bed for IGV visualization
blastn -num_threads 30 -task "blastn" -outfmt 6 -db $REF -query  $DIR/Mo17_${Chr}_${START}-${END}.fa -out $DIR/BLAST_Mo17_${Chr}_${START}-${END}_v_CI66_K10L2.out

grep "K10L2" $DIR/BLAST_Mo17_${Chr}_${START}-${END}_v_CI66_K10L2.out > $DIR/BLAST_Mo17_${Chr}_${START}-${END}_v_CI66_K10L2_K10L2.out
