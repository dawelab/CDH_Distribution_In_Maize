#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH -J TasselDemultiplexBWAaln_AllData_K10L2Contigs
#SBATCH --output TasselDemultiplexBWAaln_AllData_K10L2Contigs.out
#SBATCH --mem=600GB
#SBATCH --time=24:00:00
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
OUT_DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_K10L2Contigs"
READ_DIR="/scratch/mjb51923/raw_reads/GBS/Swarts_AllControls_LengthFiltFakeBarcodedRomeroNavarro_Romay"
KEY="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/SwartsAllControlsLengthFiltRomeroNavarroRomay_Key_SGEOnly.txt"
REF="/scratch/mjb51923/CI66_Assembly/out/CI66_rq99.asm.bp.p_ctg.gfa.fasta"
NAME="SGEOnly_K10L2"

#This defines the ref name
REF_NAME="K10L2Contigs"

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
echo "#####Converting Alignment to Data Base"xa
run_pipeline.pl -Xms200G -Xmx600G -fork1 -SAMToGBSdbPlugin -i $OUT_DIR/"BWAaln_"$NAME"_v_"$REF_NAME".sam" -db $OUT_DIR/$NAME".db" -aProp 0.0 -aLen 0 -endPlugin -runfork1

#Call SNPs
run_pipeline.pl -Xms300G -Xmx600G -fork1 -DiscoverySNPCallerPluginV2 -db $OUT_DIR/$NAME".db" -ref $REF -endPlugin -runfork1

#Asses SNP quality
run_pipeline.pl -Xms300G -Xmx600G -fork1 -SNPQualityProfilerPlugin -db $OUT_DIR/$NAME".db" -statFile $OUT_DIR/$NAME"_SNP_Qual.txt" -endPlugin -runfork1

###ProductionSNPCallerPluginV2 goes next 
run_pipeline.pl -Xms300G -Xmx600G -fork1 -ProductionSNPCallerPluginV2 -db $OUT_DIR/$NAME".db" -i $READ_DIR -k $KEY -o $OUT_DIR/$NAME".vcf" -e ApeKI -endPlugin -runfork1
