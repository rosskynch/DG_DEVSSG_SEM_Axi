prob_num=31
###############################################################################################################################################################
# Temporal convergence
###############################################################################################################################################################
we=0.3
# REPLACED due to issue with CFL ??
#K=26
#N=16
K=20
N=8

echo '#!/bin/bash' > oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs
echo '#PBS' -q SMP_queue >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs
echo '#PBS' -N oldb${prob_num} >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs
echo '#PBS' -o oldb${prob_num}_tempconv_we${we}_K${K}_N${N}-out.txt >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs
echo '#PBS' -e oldb${prob_num}_tempconv_we${we}_K${K}_N${N}-err.txt >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs
echo '#PBS' -l select=1:ncpus=16:mpiprocs=1 >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs
echo '#PBS' -l place=scatter:excl >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs
echo '#PBS' -l walltime=48:00:00 >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs
echo >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs
echo '#PROJECT=PR58' >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs
echo 'module load compiler/intel-13' >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs
echo 'cd $PBS_O_WORKDIR' >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs
echo >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs

cat ../../scripts/oldb_fixed/oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.sh >> oldb${prob_num}_tempconv_we${we}_K${K}_N${N}.pbs

###############################################################################################################################################################
# Spatial Convergence.
###############################################################################################################################################################
t=0.0001
for K in 7 10 20 26 31
do
echo '#!/bin/bash' > oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs
echo '#PBS' -q SMP_queue >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs
echo '#PBS' -N oldb${prob_num} >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs
echo '#PBS' -o oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}-out.txt >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs
echo '#PBS' -e oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}-err.txt >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs
echo '#PBS' -l select=1:ncpus=16:mpiprocs=1 >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs
echo '#PBS' -l place=scatter:excl >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs
echo '#PBS' -l walltime=48:00:00 >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs
echo >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs
echo '#PROJECT=PR58' >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs
echo 'module load compiler/intel-13' >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs
echo 'cd $PBS_O_WORKDIR' >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs
echo >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs

cat ../../scripts/oldb_fixed/oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.sh >> oldb${prob_num}_spatconv_we${we}_dt${t}_K${K}.pbs
done

###############################################################################################################################################################
# Looking at drag for each We, concentrate on higher K & N, at fixed dt=10^-4.
###############################################################################################################################################################
t=0.0001
for we in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
for K in 26 31
do
for N in 14 15 16 17 18
do
echo '#!/bin/bash' > oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs
echo '#PBS' -q SMP_queue >> oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs
echo '#PBS' -N oldb${prob_num} >> oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs
echo '#PBS' -o oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}-out.txt >> oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs
echo '#PBS' -e oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}-err.txt >> oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs
echo '#PBS' -l select=1:ncpus=16:mpiprocs=1 >> oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs
echo '#PBS' -l place=scatter:excl >> oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs
echo '#PBS' -l walltime=48:00:00 >> oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs
echo >> oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs
echo '#PROJECT=PR58' >> oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs
echo 'module load compiler/intel-13' >> oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs
echo 'cd $PBS_O_WORKDIR' >> oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs
echo >> oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs

cat ../../scripts/oldb_fixed/oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.sh >> oldb${prob_num}_we${we}_dt${t}_K${K}_N${N}.pbs

done
done
done
