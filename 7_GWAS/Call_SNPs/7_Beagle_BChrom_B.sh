#Define variables
DIR=""
VCF="GWAS_Mo17_Bchr.B.landrace.filt"
CDH="B"

cd $DIR

#This performs imputation for beagle.
java -Xmx90g -jar /apps/eb/Beagle/5.4.22Jul22.46e-Java-11/beagle.jar \
gt=$DIR/$VCF".vcf.gz" \
out=$DIR/$VCF".imputed" \
ne=1000 \
nthreads=24

tabix -p vcf $DIR/$VCF".imputed.vcf.gz"

