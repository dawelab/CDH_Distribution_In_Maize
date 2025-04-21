#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH -J Tassel_SNPs_AllData_B73Ab10HiFiBChromv2_v2
#SBATCH --output Tassel_SNPs_AllData_B73Ab10HiFiBChromv2_v2.out3
#SBATCH --mem=900GB
#SBATCH --time=24:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=END

#Load modules
module load TASSEL/5.2.44-Java-1.8.0_241

#Define Variables
OUT_DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"
READ_DIR="/scratch/mjb51923/raw_reads/GBS/Swarts_AllControls_LengthFiltFakeBarcodedRomeroNavarro_Romay"
KEY="/scratch/mjb51923/Ab10_Global_Survey/out/SwartsAllControlsLengthFiltRomeroNavarroRomay_Key.txt"
REF="/scratch/mjb51923/ref_genomes/Ab10_HiFi_v2_corrected_BChrom.noscaf.fa"
NAME="AllData_v5"

#Call SNPs
#run_pipeline.pl -Xms300G -Xmx850G -fork1 -DiscoverySNPCallerPluginV2 -db $OUT_DIR/$NAME".db" -ref $REF -endPlugin -runfork1

#Asses SNP quality
#run_pipeline.pl -Xms300G -Xmx850G -fork1 -SNPQualityProfilerPlugin -db $OUT_DIR/$NAME".db" -statFile $OUT_DIR/$NAME"_SNP_Qual.txt" -endPlugin -runfork1

###ProductionSNPCallerPluginV2 goes next 
run_pipeline.pl -Xms300G -Xmx850G -fork1 -ProductionSNPCallerPluginV2 -db $OUT_DIR/$NAME".db" -i $READ_DIR -k $KEY -o $OUT_DIR/$NAME".vcf" -e ApeKI -endPlugin -runfork1
