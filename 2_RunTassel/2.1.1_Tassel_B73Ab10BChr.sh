#Load modules
module load TASSEL/5.2.44-Java-1.8.0_241
module load BWA/0.7.17-GCCcore-11.3.0
module load SAMtools/0.1.20-GCC-11.3.0
module load BEDTools/2.29.2-GCC-8.3.0

#Define Variables
OUT_DIR=""
READ_DIR="/path/to/all/gbs/reads/properlly/formatted"
#This file is available under 1.5
KEY="SwartsAllControlsLengthFiltRomeroNavarroRomay_Key.txt"
REF="B73_Ab10_v2_BChrom.fa"
NAME="AllData_v5"
REF_NAME="B73-Ab10_BChrom"

#Demultiplex the Additional Control reads
echo "#####Demultiplexing Fastq Files"
run_pipeline.pl -Xms200G -Xmx600G -fork1 -GBSSeqToTagDBPlugin -e ApeKI -i $READ_DIR -db $OUT_DIR/$NAME".db" -k $KEY -kmerLength 64 -minKmerL 20 -mxKmerNum 100000000 -mnQS 0 -endPlugin -runfork1

#Index the references
#echo "#####Indexing the Reference"
bwa index -a bwtsw $REF

#Convert the database created in the previous step to a fastq file that can be interpreted by an aligner 
echo "#####Converting to Fastq File"
run_pipeline.pl -Xms200G -Xmx600G -fork1 -TagExportToFastqPlugin -db  $OUT_DIR/$NAME".db" -o $OUT_DIR/$NAME".fastq" -c 1 -endPlugin -runfork1
#Index the references
echo "#####Indexing the Reference"
bwa index -a bwtsw $REF

#Align the GBS .fastq file using BWA.
echo "#####Aligning the Reads"
bwa aln -t 30 $REF $OUT_DIR/$NAME".fastq" > $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".sai"
bwa samse $REF $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".sai" $OUT_DIR/$NAME".fastq" > $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".sam"

#Align the Swarts Takuno GBS .fastq file using BWA.
echo "#####Aligning the Reads"
bwa aln -t 30 $REF $OUT_DIR/$NAME".fastq" > $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".sai"
bwa samse $REF $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".sai" $OUT_DIR/$NAME".fastq" > $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".sam"

#Convert the SAM file produced in the step before to a GBS database. The original database provided was the one that I made earlier with the GBSSeqToTagDBPlugin
echo "#####Converting Alignment to Data Base"
run_pipeline.pl -Xms200G -Xmx600G -fork1 -SAMToGBSdbPlugin -i $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".sam" -db $OUT_DIR/$NAME".db" -aProp 0.0 -aLen 0 -endPlugin -runfork1

#Get the Distribution of Tags Present Across Taxa
echo "#####Getting Tag Taxa Distribution"
run_pipeline.pl -Xms200G -Xmx600G -fork1 -GetTagTaxaDistFromDBPlugin -db $OUT_DIR/$NAME".db" -o $OUT_DIR/"Tassel_TagTaxaDist_"$NAME"_v_"$REF_NAME".txt" -endPlugin -runfork1

#Export the file as a bed file 
echo "#####Creating a Bed File"
samtools view -Sb $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".sam" > $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".bam"
samtools sort $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".bam" $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".s"
bedtools bamtobed -i $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".s.bam" > $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".s.bed"
