# Script to create pbs files for all runs of Oldroyd B model problem numbers 33 and 34 (waters&king start-up solutions)
# with no devss on a single processor on any available node in the SMP queue
# Creates over 1500 pbs files. 

# fixed parameters
# To change the queue or walltime, edit below
prob_num=33
re=1.0
beta=0.1111111111111111111111111111111111111111
place_type=scatter
walltime=48:00:00
reqd_queue=SMP_queue
nodes=1
numcpus=1
mpicpus=1


# Calculate how many jobs will be required for the batch.
counter=0
for we in 1.0 10.0 100.0
do
for N in 4 8 12 16
do
for K in K1 K2_Kx1Ky2 K2_Kx2Ky1 K4_Kx2Ky2
do
for L in _L8 _L16 _L32 _L64
do
for t in 0.1 0.01 0.001 0.0001
do
counter=$(( counter + 1 ))
done
done
done
done
done

total_jobs=${counter}





for beta_s in 0.0 0.5 1.0
do

if [ "${prob_num}" = "33" ]; then
	in_folder=2dchannel
elif [ "${prob_num}" = "34" ]; then
	in_folder=3dpipe
fi

for we in 1.0 10.0 100.0
do

for N in 4 8 12 16
do

for K in K1 K2_Kx1Ky2 K2_Kx2Ky1 K4_Kx2Ky2
do

for L in _L8 _L16 _L32 _L64
do

infile=${inpath}${K}${L}.dat

for t in 0.1 0.01 0.001 0.0001
do

outpath=${CODEPATH}/output/oldb_fixed/oldb${prob_num}/we${we}/
outfile=${outpath}oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}


echo '#!/bin/bash' > oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo '#PBS' -q ${reqd_queue} >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo '#PBS' -N oldb${prob_num} >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo '#PBS' -o oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}-out.txt >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo '#PBS' -e oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}-err.txt >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo '#PBS' -l select=${nodes}:ncpus=${numcpus}:mpiprocs=${mpicpus} >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo '#PBS' -l place=${place_type} >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo '#PBS' -l walltime=${walltime}>> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo '#PROJECT=PR58' >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo 'module load compiler/intel-13' >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
#echo 'cd $PBS_O_WORKDIR' >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs


echo 'CODEPATH=$HOME/FINAL/' >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo 'JOBID=$PBS_JOBID' >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
#if needed:
echo 'INDEX=$PBS_ARRAY_INDEX'  >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs

echo 'WDPATH=/scratch/sacrmk/run-directory/${JOBID}'  >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo if [ ! -d $WDPATH ]; then >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo    mkdir -p $WDPATH  >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo fi >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs

#
# copy input files from home onto specified scratch directory
#
echo cp -f -R -L ${CODEPATH}/input/* ${WDPATH}/input/ >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo cp -f ${CODEPATH}/*.mod ${WDPATH}/ >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo cp -f ${CODEPATH}/*.o ${WDPATH}/ >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs
echo cp -f ${CODEPATH}/vesemd ${WDPATH}/ >> oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.pbs


# Change to scratch directory
echo cd ${WDPATH}


echo ./vesemd -input ${infile} -output output -prob ${prob_num} -N ${N} -alphaz 1.0 -deltat ${t} -re ${re} -we ${we} -beta ${beta} '>' screenoutput.out

echo Finished We=${we} N=${N} ${K}${L} Deltat=${t}

# Copy output into home directory
echo cp output_tec.dat ${outfile}_tec.dat 
echo cp output_tec_fine.dat ${outfile}_tec_fine.dat
echo cp screenoutput.out ${outfile}.out

echo rm -rf $WDPATH

done
done
done
done
done
done
