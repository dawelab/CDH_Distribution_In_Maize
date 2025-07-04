#Load modules
module load TASSEL/5.2.44-Java-1.8.0_241
module load BWA/0.7.17-GCCcore-11.3.0
module load SAMtools/0.1.20-GCC-11.3.0
module load BEDTools/2.29.2-GCC-8.3.0

#Define Variables
OUT_DIR=""
READ_DIR="Swarts_AllControls_LengthFiltFakeBarcodedRomeroNavarro_Romay"
KEY="SwartsAllControlsLengthFiltRomeroNavarroRomay_Key_GWASSNPs.txt"
REF="Zm-Mo17-REFERENCE-CAU-2.0_Ab10.fa"
NAME="GWAS_Mo17_Ab10"

#This v2 refers to the corrected Hifi genome
REF_NAME="Mo17_Ab10"

#Demultiplex the Additional Control reads
echo "#####Demultiplexing Fastq Files"
run_pipeline.pl -Xms200G -Xmx850G -fork1 -GBSSeqToTagDBPlugin -e ApeKI -i $READ_DIR -db $OUT_DIR/$NAME".db" -k $KEY -kmerLength 64 -minKmerL 20 -mxKmerNum 100000000 -mnQS 0 -endPlugin -runfork1

#Index the references
echo "#####Indexing the Reference"
bwa index -a bwtsw $REF

#Convert the database created in the previous step to a fastq file that can be interpreted by an aligner 
echo "#####Converting to Fastq File"
run_pipeline.pl -Xms200G -Xmx850G -fork1 -TagExportToFastqPlugin -db  $OUT_DIR/$NAME".db" -o $OUT_DIR/$NAME".fastq" -c 1 -endPlugin -runfork1
#Index the references
echo "#####Indexing the Reference"
bwa index -a bwtsw $REF

#Align the GBS .fastq file using BWA.
echo "#####Aligning the Reads"
bwa aln -t 30 $REF $OUT_DIR/$NAME".fastq" > $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".sai"
bwa samse $REF $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".sai" $OUT_DIR/$NAME".fastq" > $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".sam"

#Filter the SAM file to mapq 20
echo "#####Filtering SAM file"
samtools view -q 20 -o $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME"_filtered.sam" -S $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".sam"

#Convert the SAM file produced in the step before to a GBS database. The original database provided was the one that I made earlier with the GBSSeqToTagDBPlugin
echo "#####Converting Alignment to Data Base"
run_pipeline.pl -Xms200G -Xmx850G -fork1 -SAMToGBSdbPlugin -i $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME"_filtered.sam" -db $OUT_DIR/$NAME".db" -aProp 0.0 -aLen 0 -endPlugin -runfork1

#Call SNPs
echo "#####Calling SNPs"
run_pipeline.pl -Xms300G -Xmx850G -fork1 -DiscoverySNPCallerPluginV2 -db $OUT_DIR/$NAME".db" -ref $REF -endPlugin -runfork1

#Asses SNP quality
echo "#####Assessing SNP quality"
run_pipeline.pl -Xms300G -Xmx850G -fork1 -SNPQualityProfilerPlugin -db $OUT_DIR/$NAME".db" -statFile $OUT_DIR/$NAME"_SNP_Qual.txt" -endPlugin -runfork1

###ProductionSNPCallerPluginV2 goes next 
echo "#####Production SNP caller"
run_pipeline.pl -Xms300G -Xmx850G -fork1 -ProductionSNPCallerPluginV2 -db $OUT_DIR/$NAME".db" -i $READ_DIR -k $KEY -o $OUT_DIR/$NAME".vcf" -e ApeKI -endPlugin -runfork1
