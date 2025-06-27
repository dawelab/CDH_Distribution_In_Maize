#!/usr/bin/bash
#SBATCH --partition=batch
#SBATCH -J GZip_asc.sh
#SBATCH --output GZip_asc.out
#SBATCH --mem=50GB
#SBATCH --time=10:00:00
#SBATCH	--nodes=1
#SBATCH	--ntasks=5
#SBATCH --mail-user=meghan.brady@uga.edu
#SBATCH --mail-type=BEGIN,END

cd /scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/ClimateData/MaxEntData

gzip wc2.1_30s_prec_11.asc
gzip wc2.1_30s_srad_03.asc
gzip wc2.1_30s_srad_04.asc
gzip wc2.1_30s_srad_12.asc
gzip wc2.1_30s_tavg_01.asc
gzip wc2.1_30s_tavg_06.asc
gzip wc2.1_30s_tavg_07.asc
gzip wc2.1_30s_tavg_08.asc
gzip wc2.1_30s_tavg_09.asc
gzip wc2.1_30s_tavg_10.asc
gzip wc2.1_30s_tavg_11.asc
gzip wc2.1_30s_tmax_02.asc
gzip wc2.1_30s_tmax_03.asc
gzip wc2.1_30s_tmax_04.asc
gzip wc2.1_30s_tmax_05.asc
gzip wc2.1_30s_tmax_12.asc
gzip wc2.1_30s_tmin_01.asc
gzip wc2.1_30s_tmin_06.asc
gzip wc2.1_30s_tmin_07.asc
gzip wc2.1_30s_tmin_08.asc
gzip wc2.1_30s_tmin_09.asc
gzip wc2.1_30s_tmin_10.asc
gzip wc2.1_30s_tmin_11.asc
gzip wc2.1_30s_vapr_02.asc
gzip wc2.1_30s_vapr_03.asc
gzip wc2.1_30s_vapr_04.asc
gzip wc2.1_30s_vapr_05.asc
gzip wc2.1_30s_vapr_12.asc
gzip wc2.1_30s_wind_03.asc
gzip wc2.1_30s_wind_04.asc
gzip wc2.1_30s_wind_12.asc