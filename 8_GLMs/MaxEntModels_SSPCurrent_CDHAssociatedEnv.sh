#This is the raw code for the Maize MaxEnt Model projecting to current global climate data

#This is how it actually works and this is max memory

ml Java/17.0.6
cd /scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/maxent
mkdir /scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/K10L2_MaxEntStand7_SSPCurrent
java -mx40000000m -jar maxent.jar nowarnings noprefixes -E "" -E K10L2 responsecurves jackknife outputformat=logistic outputdirectory=/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/K10L2_MaxEntStand7_SSPCurrent projectionlayers=/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/MaxEnt_EnvData_Subset/projection  samplesfile=/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/SGE_GPS_MaxEnt_Format.csv environmentallayers=/scratch/mjb51923/Ab10_Global_Survey/out/WorldClim2/MaxEnt_EnvData_Subset/CDH_Association_EnvVar/K10L2 randomseed noaskoverwrite randomtestpoints=25 betamultiplier=3.0 replicates=10 replicatetype=subsample threshold writeplotdata noautofeature nooutputgrids threads=24 defaultprevalence=0.12 -t sq5.stand -N wc2.1_30s_vapr_avg.stand autorun
