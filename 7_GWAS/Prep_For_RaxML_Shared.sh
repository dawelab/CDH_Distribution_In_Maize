#!/usr/bin/bash
#SBATCH --partition=batch 
#SBATCH -J Prep_For_RaxML_Shared
#SBATCH --output Prep_For_RaxML_Shared.out
#SBATCH --mem=150GB
#SBATCH --time=72:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=30
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load BCFtools/1.15.1-GCC-11.3.0
module load Clustal-Omega/1.2.4-GCC-11.2.0
module load SAMtools/1.18-GCC-12.3.0
module load RAxML-NG/1.2.2-GCC-12.2.0

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"

#The output from vcfF is gzipped. I have to unzip it
gunzip $DIR/AllData_v5.Ab10Shared.ControlsST.filt.vcf.gz
#I need to recompress it using the bgzip command for compatibility with bcftools
bgzip $DIR/AllData_v5.Ab10Shared.ControlsST.filt.vcf
#I index the file 
bcftools index $DIR/AllData_v5.Ab10Shared.ControlsST.filt.vcf.gz

while read line; do
#This creates a consensus sequence for one sample at a time
bcftools consensus -f /scratch/mjb51923/ref_genomes/Ab10_HiFi_v2_corrected_BChrom.noscaf.chr10.fa -s $line -o $DIR/Consensus/$line.tmp.fa $DIR/AllData_v5.Ab10Shared.ControlsST.filt.vcf.gz
#This selects only the Ab10 haplotype or one shared region
samtools faidx $DIR/Consensus/$line.tmp.fa 10:141112544-195055488 > $DIR/Consensus/$line"_Ab10Hap.fa"
samtools faidx $DIR/Consensus/$line.tmp.fa 10:141112544-142472000 > $DIR/Consensus/$line"_Shared1.fa"
samtools faidx $DIR/Consensus/$line.tmp.fa 10:152050000-156880000 > $DIR/Consensus/$line"_Shared2.fa"
samtools faidx $DIR/Consensus/$line.tmp.fa 10:158250000-170000000 > $DIR/Consensus/$line"_Shared3.fa"

#This makes the name of the sequence the same as the sample name instead of just 10
sed -i "s/>10:141112544-195055488/>$line/g" $DIR/Consensus/$line"_Ab10Hap.fa"
sed -i "s/>10:141112544-142472000/>$line/g" $DIR/Consensus/$line"_Shared1.fa"
sed -i "s/>10:152050000-156880000/>$line/g" $DIR/Consensus/$line"_Shared2.fa"
sed -i "s/>10:158250000-170000000/>$line/g" $DIR/Consensus/$line"_Shared3.fa"
done < $DIR/HigginsN10.txt

#This removes temporary files
rm *.tmp.fa

#This concatenates all the fasta files into a single multi sequence file
cat $DIR/Consensus/*_Ab10Hap.fa > $DIR/Consensus/HigginsN10_Ab10Hap.fa
cat $DIR/Consensus/*_Shared1.fa > $DIR/Consensus/HigginsN10_Shared1.fa
cat $DIR/Consensus/*_Shared2.fa > $DIR/Consensus/HigginsN10_Shared2.fa
cat $DIR/Consensus/*_Shared3.fa > $DIR/Consensus/HigginsN10_Shared3.fa

#This command generates a ML tree using the GTR+G model in RAxML
raxml-ng --check --all --msa $DIR/HigginsN10_Shared3.aln.fasta --model "GTR+G" --seed 10 --prefix HigginsN10_Shared3 --threads 10
raxml-ng --check --all --msa $DIR/HigginsN10_Shared2.aln.fasta --model "GTR+G" --seed 10 --prefix HigginsN10_Shared2 --threads 10
raxml-ng --check --all --msa $DIR/HigginsN10_Shared1.aln.fasta --model "GTR+G" --seed 10 --prefix HigginsN10_Shared1 --threads 10
raxml-ng --check --all --msa $DIR/HigginsN10_Ab10Hap.aln.fasta --model "GTR+G" --seed 10 --prefix HigginsN10_Ab10Hap --threads 10