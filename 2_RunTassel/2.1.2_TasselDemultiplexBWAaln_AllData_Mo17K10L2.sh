#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH -J TasselDemultiplexBWAaln_AllData_Mo17K10L2
#SBATCH --output TasselDemultiplexBWAaln_AllData_Mo17K10L2.out
#SBATCH --mem=600GB
#SBATCH --time=48:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=30
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=END

#Load modules
module load TASSEL/5.2.44-Java-1.8.0_241
module load BWA/0.7.17-GCCcore-11.3.0
module load SAMtools/0.1.20-GCC-11.3.0
module load BEDTools/2.29.2-GCC-8.3.0

#Define Variables
OUT_DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_Mo17K10L2"
READ_DIR="/scratch/mjb51923/raw_reads/GBS/Swarts_AllControls_LengthFiltFakeBarcodedRomeroNavarro_Romay"
KEY="/scratch/mjb51923/Ab10_Global_Survey/out/SwartsAllControlsLengthFiltRomeroNavarroRomay_Key.txt"
REF="/scratch/mjb51923/ref_genomes/Zm-Mo17-REFERENCE-CAU-2.0.noscaf.K10L2.fa"
NAME="AllData_v6"
REF_NAME="Mo17K10L2"

#Demultiplex the Additional Control reads
echo "#####Demultiplexing Fastq Files"
run_pipeline.pl -Xms200G -Xmx600G -fork1 -GBSSeqToTagDBPlugin -e ApeKI -i $READ_DIR -db $OUT_DIR/$NAME".db" -k $KEY -kmerLength 64 -minKmerL 20 -mxKmerNum 100000000 -mnQS 0 -endPlugin -runfork1

#Index the references
#echo "#####Indexing the Reference"
#bwa index -a bwtsw $REF

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
