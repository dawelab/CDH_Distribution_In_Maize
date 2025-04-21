#!/usr/bin/bash
#SBATCH --partition=batch 
#SBATCH -J Download_FAO_SoilData
#SBATCH --output Download_FAO_SoilData.out
#SBATCH --mem=5GB
#SBATCH --time=10:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=1
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

cd /scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/SoilData

wget https://www.fao.org/fileadmin/user_upload/soils/docs/HWSD/Soil_Quality_data/sq1.asc
wget https://www.fao.org/fileadmin/user_upload/soils/docs/HWSD/Soil_Quality_data/sq2.asc
wget https://www.fao.org/fileadmin/user_upload/soils/docs/HWSD/Soil_Quality_data/sq3.asc
wget https://www.fao.org/fileadmin/user_upload/soils/docs/HWSD/Soil_Quality_data/sq3.asc
wget https://www.fao.org/fileadmin/user_upload/soils/docs/HWSD/Soil_Quality_data/sq4.asc
wget https://www.fao.org/fileadmin/user_upload/soils/docs/HWSD/Soil_Quality_data/sq5.asc
wget https://www.fao.org/fileadmin/user_upload/soils/docs/HWSD/Soil_Quality_data/sq6.asc
wget https://www.fao.org/fileadmin/user_upload/soils/docs/HWSD/Soil_Quality_data/sq7.asc