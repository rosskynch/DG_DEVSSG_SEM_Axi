for prob_num in 33 34
do
for N in 4 8 12 16
do
for we in 1.0 10.0 100.0
do
for L in _L8 _L16 _L32 _L64
do
echo '#!/bin/bash' > oldb${prob_num}_we${we}${L}_N${N}.pbs
echo '#PBS' -q workq >> oldb${prob_num}_we${we}${L}_N${N}.pbs
echo '#PBS' -N oldb${prob_num} >> oldb${prob_num}_we${we}${L}_N${N}.pbs
echo '#PBS' -o oldb${prob_num}_we${we}${L}_N${N}-out.txt >> oldb${prob_num}_we${we}${L}_N${N}.pbs
echo '#PBS' -e oldb${prob_num}_we${we}${L}_N${N}-err.txt >> oldb${prob_num}_we${we}${L}_N${N}.pbs
echo '#PBS' -l select=1:ncpus=16:mpiprocs=16 >> oldb${prob_num}_we${we}${L}_N${N}.pbs
echo '#PBS' -l place=scatter:excl >> oldb${prob_num}_we${we}${L}_N${N}.pbs
echo '#PBS' -l walltime=72:00:00 >> oldb${prob_num}_we${we}${L}_N${N}.pbs
echo >> oldb${prob_num}_we${we}${L}_N${N}.pbs
echo '#PROJECT=PR58' >> oldb${prob_num}_we${we}${L}_N${N}.pbs
echo 'module load compiler/intel-13' >> oldb${prob_num}_we${we}${L}_N${N}.pbs
echo 'cd $PBS_O_WORKDIR' >> oldb${prob_num}_we${we}${L}_N${N}.pbs
echo >> oldb${prob_num}_we${we}${L}_N${N}.pbs

cat ../../scripts/oldb_fixed/oldb${prob_num}_we${we}${L}_N${N}.sh >> oldb${prob_num}_we${we}${L}_N${N}.pbs
done
done
done
done
