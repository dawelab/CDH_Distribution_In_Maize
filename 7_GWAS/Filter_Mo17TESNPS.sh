#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Filter_Mo17TESNPS
#SBATCH --output Filter_Mo17TESNPS.out
#SBATCH --mem=100GB
#SBATCH --time=24:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=30
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

module load BEDTools/2.31.0-GCC-12.3.0
module load BCFtools/1.18-GCC-12.3.0

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/"
VCF="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/SGEOnly_Mo17.chr.filt.vcf.gz"
TE="/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/Zm-Mo17-REFERENCE-CAU-2.0.TE.gff3.gz"
OUT="SGEOnly_Mo17.chr.filtTE.vcf"

bedtools intersect -v -a $VCF -b $TE > $DIR/$OUT

bgzip $DIR/$OUT

DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/"
VCF="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2/SGEOnly_Mo17.chr.filt.landrace.vcf.gz"
TE="/scratch/mjb51923/Ab10_Global_Survey/out/OrthoFinder/Zm-Mo17-REFERENCE-CAU-2.0.TE.gff3.gz"
OUT="SGEOnly_Mo17.chr.landrace.filtTE.vcf"

bedtools intersect -v -a $VCF -b $TE > $DIR/$OUT

bgzip $DIR/$OUT