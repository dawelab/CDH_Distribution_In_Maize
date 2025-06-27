#Load modules
module load Beagle/5.4.22Jul22.46e-Java-11
module load BCFtools/1.15.1-GCC-11.3.0

#Define variables
DIR=""
VCF="GWAS_Mo17_Ab10.Ab10.landrace.filt"
CDH="Ab10"

cd $DIR

#This performs imputation for beagle. 
java -Xmx90g -jar /apps/eb/Beagle/5.4.22Jul22.46e-Java-11/beagle.jar \
gt=$DIR/$VCF".vcf.gz" \
out=$DIR/$VCF".imputed" \
ne=1000 \
nthreads=24

tabix -p vcf $DIR/$VCF".imputed.vcf.gz"


