# Fixed:
prob_num=31
re=0.0
beta=0.59
beta_s=1.0
inpath=input/fixedcylinder/fixedcylindermeshK 

# Fixed for convergence checks:
we=0.3

###############################################################################################################################################################
# Temporal convergence
# Switched to run 4 at once.
# The highest deltat should take the longest to run, so we shouldn't lose any time on the lower values.
###############################################################################################################################################################
# Note: @K=26 N =16 - dt=0.01 didn't work.. CFL condition ?? - try 0.005 instead?
# Instead, trying to lower K/N and run more timesteps - may convergence or a "wrong value" ??
K=20
N=8
outpath=output/oldb_devss1.0/oldb${prob_num}/tempconv/ 

echo '#Temporal convergence script for prob'${prob_num} > oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.sh
echo 'export OMP_NUM_THREADS=1' >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.sh
echo 'cd ../../' >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.sh
echo >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.sh

for t in 0.01 0.001 0.0001 0.00001 
do

outfile_screen=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.out
outfile=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}
infile=${inpath}${K}.dat

echo ./vesemd -input ${infile} -output ${outfile} -prob ${prob_num} -N ${N} -alphaz 1.0 -deltat ${t} -re ${re} -we ${we} -beta ${beta} -beta_s ${beta_s} '>' ${outfile_screen}' && echo Finished 'We=${we} N=${N} K=${K} Deltat=${t} '&' >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.sh

done
echo 'wait' >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.sh

###############################################################################################################################################################
# Spatial Convergence.
# Switched to run 4 per mesh.
# The highest N should take the longest to run, so we shouldn't lose any time on the lower values.
###############################################################################################################################################################
t=0.0001
outpath=output/oldb_devss1.0/oldb${prob_num}/spatconv/

for K in 7 10 20 26 31
do
echo '#Spatial convergence script for prob'${prob_num} > oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.sh
echo 'export OMP_NUM_THREADS=1' >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.sh
echo 'cd ../../' >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.sh
echo >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.sh

for N in 4 8 12 16
do

outfile_screen=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.out
outfile=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}
infile=${inpath}${K}.dat

echo ./vesemd -input ${infile} -output ${outfile} -prob ${prob_num} -N ${N} -alphaz 1.0 -deltat ${t} -re ${re} -we ${we} -beta ${beta} -beta_s ${beta_s} '>' ${outfile_screen}' && echo Finished 'We=${we} N=${N} K=${K} Deltat=${t} '&' >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.sh


done
echo 'wait' >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.sh

done

###############################################################################################################################################################
# Looking at drag for each We, concentrate on higher K & N, at fixed dt=10^-4.
# Changed to run 4 instances at once per mesh per weissenberg number.
# The highest N should take the longest to run, so we shouldn't lose any time on the lower values.
###############################################################################################################################################################
t=0.0001

for we in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
outpath=output/oldb_devss1.0/oldb${prob_num}/we${we}/
for K in 26 31
do

if [ "$K" = "26" ]; then

echo '# Script for prob'${prob_num} we=${we} dt=${t} K=${K} > oldb${prob_num}_we${we}_dt${t}_K${K}.sh
echo 'export OMP_NUM_THREADS=1' >> oldb${prob_num}_we${we}_dt${t}_K${K}.sh
echo 'cd ../../' >> oldb${prob_num}_we${we}_dt${t}_K${K}.sh
echo >> oldb${prob_num}_we${we}_dt${t}_K${K}.sh

for N in 15 16 17 18
do

outfile_screen=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.out
outfile=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}
infile=${inpath}${K}.dat

echo ./vesemd -input ${infile} -output ${outfile} -prob ${prob_num} -N ${N} -alphaz 1.0 -deltat ${t} -re ${re} -we ${we} -beta ${beta} -beta_s ${beta_s} '>' ${outfile_screen}' && echo Finished 'We=${we} N=${N} K=${K} Deltat=${t} '&' >> oldb${prob_num}_we${we}_dt${t}_K${K}.sh

done
echo 'wait' >> oldb${prob_num}_we${we}_dt${t}_K${K}.sh

elif [ "$K" = "31" ]; then

echo '# Script for prob'${prob_num} we=${we} dt=${t} K=${K} > oldb${prob_num}_we${we}_dt${t}_K${K}.sh
echo 'export OMP_NUM_THREADS=1' >> oldb${prob_num}_we${we}_dt${t}_K${K}.sh
echo 'cd ../../' >> oldb${prob_num}_we${we}_dt${t}_K${K}.sh
echo >> oldb${prob_num}_we${we}_dt${t}_K${K}.sh

for N in 13 14 15 16
do

outfile_screen=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.out
outfile=${outpath}oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}
infile=${inpath}${K}.dat

echo ./vesemd -input ${infile} -output ${outfile} -prob ${prob_num} -N ${N} -alphaz 1.0 -deltat ${t} -re ${re} -we ${we} -beta ${beta} -beta_s ${beta_s} '>' ${outfile_screen}' && echo Finished 'We=${we} N=${N} K=${K} Deltat=${t} '&' >> oldb${prob_num}_we${we}_dt${t}_K${K}.sh

done
echo 'wait' >> oldb${prob_num}_we${we}_dt${t}_K${K}.sh

fi
done
done
###############################################################################################################################################################



