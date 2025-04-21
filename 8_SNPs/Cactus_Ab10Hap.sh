#!/usr/bin/bash
#SBATCH --partition=highmem_p 
#SBATCH -J Cactus_Ab10Hap
#SBATCH --output Cactus_Ab10Hap.out
#SBATCH --mem=600GB
#SBATCH --time=72:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=30
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load Cactus/2.6.7-GCCcore-11.3.0-Python-3.10.4
module load SAMtools/1.18-GCC-12.3.0

module load numpy/1.9.2-intel-2021b-Python-2.7.18
module load Python/3.11.5-GCCcore-13.2.0

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"

#This generates the reference sequence
samtools faidx /scratch/mjb51923/ref_genomes/Ab10_HiFi_v2_corrected_BChrom.noscaf.chr10.fa 10:141112544-195055488 > $DIR/Ab10_HiFi_v2_corrected_BChrom.noscaf.Ab10Hap.fa

#This changes the sequence names to Ab10 instead of the sample name within each file
while read line; do
cp $DIR/Consensus/$line"_Ab10Hap.fa" $DIR/Consensus/$line"_Ab10Hap.Cactus.fa"
sed -i "s/>$line/>Ab10/g" $DIR/Consensus/$line"_Ab10Hap.Cactus.fa"
done < $DIR/HigginsN10.txt

#Align the sequences
quay.io/comparative-genomics-toolkit/cactus:v2.9.3
cactus-pangenome $DIR/Cactus_Ab10Hap.txt --outDir $DIR/Cactus_Ab10Hap --outName Ab10_Hap --reference $DIR/Ab10_HiFi_v2_corrected_BChrom.noscaf.Ab10Hap.fa --maxThreads 30
