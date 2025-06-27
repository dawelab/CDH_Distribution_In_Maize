#Load the module
module load BCFtools/1.15.1-GCC-11.3.0

#Define the Variables
DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"
VCF="GWAS_Mo17_Bchr"
CDH="B"

#This zips the vcf file
bgzip $DIR/$VCF.vcf

#Index the new vcf file
bcftools index $DIR/$VCF".vcf.gz"

#Subset the SNPs to only include those in the 1-10 range
bcftools view  --threads 24 -r 1,2,3,4,5,6,7,8,9,10 -o $DIR/$VCF".chr.vcf.gz" $DIR/$VCF.vcf.gz

#Index the new vcf file
bcftools index $DIR/$VCF".chr.vcf.gz"

#I need to define Quality score in the header
This isolates the header
bcftools view -h $DIR/${VCF}.chr.vcf.gz > $DIR/${CDH}_hdr.txt

#I manually added the Quality score to the header

#This adds the new header
bcftools reheader -h $DIR/${CDH}_hdr.txt $DIR/${VCF}.chr.vcf.gz > $DIR/${VCF}.chrfix.vcf.gz

#Subset to select only samples that are landraces
bcftools view --threads 24 -S $DIR/Landrace_Only.txt -o $DIR/${VCF}.chr.landrace.vcf.gz $DIR/${VCF}.chrfix.vcf.gz

#Index the new vcf files
bcftools index $DIR/${VCF}.chr.landrace.vcf.gz

#Subset to select only alleles with at least a 5% minor allele frequency, a depth of 3 to 20, and a genotype quality of 60, and have sample missingness of less than 25%
bcftools filter --threads 24 -i 'FORMAT/DP > 3  && FORMAT/DP < 20 && MAF[0] > 0.05 && GQ > 60 && F_MISSING<0.25' -o $DIR/${VCF}.chr.landrace.filt.vcf.gz $DIR/${VCF}.chr.landrace.vcf.gz

#Index the new vcf files
bcftools index $DIR/${VCF}.chr.landrace.filt.vcf.gz

######################################################
#This section repeats the subsetting for the CDH as a control
bcftools view  --threads 24 -r AB10 -o $DIR/$VCF'.'$CDH'.vcf.gz' $DIR/$VCF.vcf.gz

#Index the new vcf file
bcftools index $DIR/$VCF'.'$CDH'.vcf.gz'

#This adds the new header
bcftools reheader -h $DIR/${CDH}_hdr.txt $DIR/$VCF'.'$CDH'.vcf.gz' > $DIR/$VCF'.'$CDH'.fix.vcf.gz'

#Subset to select only samples that are landraces
bcftools view --threads 24 -S $DIR/Landrace_Only.txt -o $DIR/$VCF'.'$CDH'.landrace.vcf.gz' $DIR/$VCF'.'$CDH'.fix.vcf.gz'

#Index the new vcf file
bcftools index $DIR/$VCF'.'$CDH'.landrace.vcf.gz'

#Subset to select only alleles with at least a 5% minor allele frequency, a depth of 3 to 20, and a genotype quality of 60, and have sample missingness of less than 25%
bcftools filter --threads 24 -i 'FORMAT/DP > 3  && FORMAT/DP < 20 && MAF[0] > 0.05 && GQ > 60 && F_MISSING<0.25' -o $DIR/$VCF'.'$CDH'.landrace.filt.vcf.gz' $DIR/$VCF'.'$CDH'.landrace.vcf.gz'

#Index the new vcf file
bcftools index $DIR/$VCF'.'$CDH'.landrace.filt.vcf.gz'
