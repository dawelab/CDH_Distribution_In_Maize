#!/usr/bin/bash
#SBATCH --partition=batch 
#SBATCH -J Download_WorldClim2_Future
#SBATCH --output Download_WorldClim2_Future.out
#SBATCH --mem=5GB
#SBATCH --time=10:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

# cd /scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/ClimateData/SSP126
# for i in "tmin" \
# "tmax" \
# "prec" \
# "bioc"
# do
# wget https://geodata.ucdavis.edu/cmip6/30s/MIROC6/ssp126/wc2.1_30s_${i}_MIROC6_ssp126_2081-2100.tif
# done

# cd /scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/ClimateData/SSP245
# for i in "tmin" \
# "tmax" \
# "prec" \
# "bioc"
# do
# wget https://geodata.ucdavis.edu/cmip6/30s/MIROC6/ssp245/wc2.1_30s_${i}_MIROC6_ssp245_2081-2100.tif
# done

# cd /scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/ClimateData/SSP370
# for i in "tmin" \
# "tmax" \
# "prec" \
# "bioc"
# do
# wget https://geodata.ucdavis.edu/cmip6/30s/MIROC6/ssp370/wc2.1_30s_${i}_MIROC6_ssp370_2081-2100.tif
# done

cd /scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/ClimateData/SSP585
for i in "tmin" \
"tmax" \
"prec" \
"bioc"
do
wget https://geodata.ucdavis.edu/cmip6/30s/MIROC6/ssp585/wc2.1_30s_${i}_MIROC6_ssp585_2081-2100.tif
done
