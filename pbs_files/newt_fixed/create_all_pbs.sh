for prob_num in 21 22 23 24
do
for re in 0.001 0.01 0.1 1.0 10.0
do
echo '#!/bin/bash' > newt${prob_num}_re${re}.pbs
echo '#PBS' -q workq >> newt${prob_num}_re${re}.pbs
echo '#PBS' -N newt${prob_num} >> newt${prob_num}_re${re}.pbs
echo '#PBS' -o newt${prob_num}_re${re}-out.txt >> newt${prob_num}_re${re}.pbs
echo '#PBS' -e newt${prob_num}_re${re}-err.txt >> newt${prob_num}_re${re}.pbs
echo '#PBS' -l select=1:ncpus=8:mpiprocs=8 >> newt${prob_num}_re${re}.pbs
echo '#PBS' -l place=scatter:excl >> newt${prob_num}_re${re}.pbs
echo '#PBS' -l walltime=72:00:00 >> newt${prob_num}_re${re}.pbs
echo >> newt${prob_num}_re${re}.pbs
echo '#PROJECT=PR58' >> newt${prob_num}_re${re}.pbs
echo 'module load compiler/intel-13' >> newt${prob_num}_re${re}.pbs
echo 'cd $PBS_O_WORKDIR' >> newt${prob_num}_re${re}.pbs
echo >> newt${prob_num}_re${re}.pbs

cat ../../scripts/newt_fixed/newt${prob_num}_re${re}.sh >> newt${prob_num}_re${re}.pbs
done
done
