for beta_name in 0.1 0.5 0.9
do

if [ "$beta_name" = "0.1" ]; then
beta=0.111111111111111111111111
elif [ "$beta_name" = "0.5" ]; then
beta=0.5
elif [ "$beta_name" = "0.9" ]; then
beta=0.888888888888888888888889
fi


# Fixed:
prob_num=32
re=0.01
inpath=input/fixedsphere/fixedspheremeshK

# Fixed for convergence checks:
we=0.3

###############################################################################################################################################################
# Temporal convergence
###############################################################################################################################################################
# Note: @K=26 N =16 - dt=0.01 didn't work.. CFL condition ?? - try 0.005 instead?
# Instead, trying to lower K/N and run more timesteps - may convergence or a "wrong value" ??
K=20
N=8
outpath=output/oldb_fixed/oldb${prob_num}/beta${beta_name}/tempconv/ 

echo '#Temporal convergence script for prob'${prob_num} > oldb${prob_num}_tempconv_beta${beta_name}_we${we}_K${K}_N${N}.sh
echo 'export OMP_NUM_THREADS=1' >> oldb${prob_num}_tempconv_beta${beta_name}_we${we}_K${K}_N${N}.sh
echo 'cd ../../' >> oldb${prob_num}_tempconv_beta${beta_name}_we${we}_K${K}_N${N}.sh
echo >> oldb${prob_num}_tempconv_beta${beta_name}_we${we}_K${K}_N${N}.sh

for t in 0.01 0.001 0.0001 0.00001
do

outfile_screen=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.out
outfile=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}
infile=${inpath}${K}.dat

echo ./vesemd -input ${infile} -output ${outfile} -prob ${prob_num} -N ${N} -alphaz 1.0 -deltat ${t} -re ${re} -we ${we} -beta ${beta} '>' ${outfile_screen} >> oldb${prob_num}_tempconv_beta${beta_name}_we${we}_K${K}_N${N}.sh

echo 'echo Finished a run' >> oldb${prob_num}_tempconv_beta${beta_name}_we${we}_K${K}_N${N}.sh

done

###############################################################################################################################################################
# Spatial Convergence.
###############################################################################################################################################################
t=0.0001
outpath=output/oldb_fixed/oldb${prob_num}/beta${beta_name}/spatconv/
for K in 7 10 20 26 31
do

echo '#Spatial convergence script for prob'${prob_num} > oldb${prob_num}_spatconv_beta${beta_name}_we${we}_dt${t}_K${K}.sh
echo 'export OMP_NUM_THREADS=1' >> oldb${prob_num}_spatconv_beta${beta_name}_we${we}_dt${t}_K${K}.sh
echo 'cd ../../' >> oldb${prob_num}_spatconv_beta${beta_name}_we${we}_dt${t}_K${K}.sh
echo >> oldb${prob_num}_spatconv_beta${beta_name}_we${we}_dt${t}_K${K}.sh

for N in 4 8 12 16
do

outfile_screen=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.out
outfile=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}
infile=${inpath}${K}.dat

echo ./vesemd -input ${infile} -output ${outfile} -prob ${prob_num} -N ${N} -alphaz 1.0 -deltat ${t} -re ${re} -we ${we} -beta ${beta} '>' ${outfile_screen} >> oldb${prob_num}_spatconv_beta${beta_name}_we${we}_dt${t}_K${K}.sh

echo 'echo Finished a run' >> oldb${prob_num}_spatconv_beta${beta_name}_we${we}_dt${t}_K${K}.sh

done
done

###############################################################################################################################################################
# Looking at drag for each We, concentrate on higher K & N, at fixed dt=10^-4.
###############################################################################################################################################################
t=0.0001

for we in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3
do
outpath=output/oldb_fixed/oldb${prob_num}/beta${beta_name}/we${we}/
for K in 26 31
do
for N in 14 15 16 17 18
do

echo '# Script for prob'${prob_num} we=${we} dt=${t} K=${K} N=${N} > oldb${prob_num}_beta${beta_name}_we${we}_dt${t}_K${K}_N${N}.sh
echo 'export OMP_NUM_THREADS=1' >> oldb${prob_num}_beta${beta_name}_we${we}_dt${t}_K${K}_N${N}.sh
echo 'cd ../../' >> oldb${prob_num}_beta${beta_name}_we${we}_dt${t}_K${K}_N${N}.sh
echo >> oldb${prob_num}_beta${beta_name}_we${we}_dt${t}_K${K}_N${N}.sh

outfile_screen=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.out
outfile=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}
infile=${inpath}${K}.dat

echo ./vesemd -input ${infile} -output ${outfile} -prob ${prob_num} -N ${N} -alphaz 1.0 -deltat ${t} -re ${re} -we ${we} -beta ${beta} '>' ${outfile_screen} >> oldb${prob_num}_beta${beta_name}_we${we}_dt${t}_K${K}_N${N}.sh

echo 'echo Finished a run' >> oldb${prob_num}_beta${beta_name}_we${we}_dt${t}_K${K}_N${N}.sh

done
done
done
###############################################################################################################################################################

done



