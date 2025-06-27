module load BCFtools/1.15.1-GCC-11.3.0

DIR=""
VCF=GWAS_Mo17_Bchr.chr.landrace.filt

for chr in {1..10}; do
#This subsets to just one chromosome
  bcftools view --threads 24 -r $chr $DIR/$VCF.vcf.gz -Oz -o $DIR/${VCF}_chr${chr}.vcf.gz
  tabix -p vcf $DIR/${VCF}_chr${chr}.vcf.gz
done
