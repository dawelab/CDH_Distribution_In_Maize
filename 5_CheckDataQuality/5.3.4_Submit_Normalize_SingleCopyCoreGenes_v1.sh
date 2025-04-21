#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Submit_Normalize_SingleCopyCoreGenes_v1
#SBATCH --output Submit_Normalize_SingleCopyCoreGenes_v1.%A-%a.out
#SBATCH --mem=200GB
#SBATCH --time=100:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=12
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --array=1-20

#This defines variables
DIR="/scratch/mjb51923/Ab10_Global_Survey/out"
DIRNAME="AlignGBS_HiFiAb10Corrected_v2"
NAME="Tassel_TagTaxaDist_AllData_v5_v_B73-Ab10HIFI_B-Chrom_v2"
IT=$SLURM_ARRAY_TASK_ID
SGE="SingleCopyCoreGenes"

#This writes out what file is working 
echo "working on file"
echo $DIRNAME
echo "normalizing this SGE"
echo $SGE

#Load the modules needed for sapelo2
module load R/4.3.1-foss-2022a

Rscript --vanilla /scratch/mjb51923/Ab10_Global_Survey/scripts/Normalize_SingleCopyCoreGenes.R $DIR/$DIRNAME $NAME $IT $SGE
