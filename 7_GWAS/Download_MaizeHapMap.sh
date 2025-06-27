#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J Download_MaizeHapMap
#SBATCH --output Download_MaizeHapMap.out
#SBATCH --mem=50GB
#SBATCH --time=24:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=24
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END


#Load the module
module load BCFtools/1.15.1-GCC-11.3.0

#Define the Variables
DIR="/scratch/mjb51923/Ab10_Global_Survey/out/AlignGBS_HiFiAb10Corrected_v2"

cd $DIR

#dowload the maize haplotype map from datacommons
wget https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/hmp321_imputed/merged_flt_c10.imputed.vcf.gz
wget https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/hmp321_imputed/merged_flt_c1.imputed.vcf.gz
wget https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/hmp321_imputed/merged_flt_c2.imputed.vcf.gz
wget https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/hmp321_imputed/merged_flt_c3.imputed.vcf.gz
wget https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/hmp321_imputed/merged_flt_c4.imputed.vcf.gz
wget https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/hmp321_imputed/merged_flt_c5.imputed.vcf.gz
wget https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/hmp321_imputed/merged_flt_c6.imputed.vcf.gz
wget https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/hmp321_imputed/merged_flt_c7.imputed.vcf.gz
wget https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/hmp321_imputed/merged_flt_c8.imputed.vcf.gz
wget https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/hmp321_imputed/merged_flt_c9.imputed.vcf.gz

bcftools concat -o $DIR/merged_flt_out.imputed.vcf.gz  $DIR/merged_flt_c*.imputed.vcf.gz 