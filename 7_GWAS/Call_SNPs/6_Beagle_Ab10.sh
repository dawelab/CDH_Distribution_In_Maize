#This defines the variables
IT=$SLURM_ARRAY_TASK_ID

#Load modules
module load Beagle/5.4.22Jul22.46e-Java-11
module load BCFtools/1.15.1-GCC-11.3.0

#Define variables
DIR=""
VCF="GWAS_Mo17_Ab10.chr.landrace.filt"
CDH="AB10"

#This defines all chromosomes
CHR=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
#This isolates a single chromosome
chr=${CHR[$IT-1]}
echo "############ working on chr"$chr

#This performs imputation for beagle. I am using an effective population size of 1000 for the domestication bottle neck
java -Xmx90g -jar /apps/eb/Beagle/5.4.22Jul22.46e-Java-11/beagle.jar \
gt=$DIR/${VCF}_chr${chr}.vcf.gz \
out=$DIR/${VCF}_chr${chr}.imputed \
ne=1000 \
nthreads=24

tabix -p vcf $DIR/${VCF}_chr${chr}.imputed.vcf.gz

