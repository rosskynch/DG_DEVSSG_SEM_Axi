# Script to create pbs batch for all runs of Oldroyd B model problem numbers 33 (waters&king start-up solutions)
# with no devss on a single processor on any available node in the SMP queue

# fixed parameters
# To change the queue or walltime, edit below
prob_num=34
re=1.0
beta=0.1111111111111111111111111111111111111111
place_type=scatter
walltime=48:00:00
reqd_queue=SMP_queue
nodes=1
numcpus=1
mpicpus=1

model=oldb
output_dir=oldb_fixed
CODEPATH=$HOME/FINAL/
job_name=${model}${prob_num}_batch

if [ "${prob_num}" = "31" ]; then
	in_folder=fixedcylinder
elif [ "${prob_num}" = "32" ]; then
	in_folder=fixedsphere
elif [ "${prob_num}" = "33" ]; then
	in_folder=2dchannel
elif [ "${prob_num}" = "34" ]; then
	in_folder=3dpipe
fi

inpath=input/${in_folder}

# Calculate how many jobs will be required for the batch.
counter=0
for we in 1.0 10.0 100.0
do
for N in 4 8 12 16
do
for K in K2_Kx1Ky2
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

pbs_file=${model}${prob_num}_batch.pbs

echo '#!/bin/bash' > ${pbs_file}
echo '#PBS' -J 1-${total_jobs} >> ${pbs_file}
echo '#PBS' -q ${reqd_queue} >> ${pbs_file}
echo '#PBS' -N ${job_name} >> ${pbs_file}
#echo '#PBS' -o ${model}${prob_num}_batch-out.txt >> ${pbs_file} # Removed
#echo '#PBS' -e ${model}${prob_num}_batch-err.txt >> ${pbs_file} # Removed | this and the above change where stderr and stdout go, but with a batch of jobs, they would overwrite each other.
echo '#PBS' -l select=${nodes}:ncpus=${numcpus}:mpiprocs=${mpicpus} >> ${pbs_file}
echo '#PBS' -l place=${place_type} >> ${pbs_file}
echo '#PBS' -l walltime=${walltime} >> ${pbs_file}
echo >> ${pbs_file}
echo '#PROJECT=PR58' >> ${pbs_file}
echo 'module load compiler/intel-13' >> ${pbs_file}
echo >> ${pbs_file}


# Shorten PBS paramters.
echo 'JOBID=$PBS_JOBID' >> ${pbs_file}
echo 'INDEX=$PBS_ARRAY_INDEX'  >> ${pbs_file}

echo >> ${pbs_file}

# Set up code and working paths.
echo 'CODEPATH=$HOME/FINAL' >> ${pbs_file}
echo 'WDPATH=/scratch/sacrmk/run-directory/'${job_name}'/${INDEX}' >> ${pbs_file}
echo 'if [ ! -d $WDPATH ]; then' >> ${pbs_file}
echo 'mkdir -p $WDPATH'  >> ${pbs_file}
echo 'fi' >> ${pbs_file}

echo >> ${pbs_file}

#
# copy code and input files from home into specified scratch directory
#
echo 'cp -f -R -L ${CODEPATH}/input/ ${WDPATH}/' >> ${pbs_file}
echo 'cp -f ${CODEPATH}/*.mod ${WDPATH}/' >> ${pbs_file}
echo 'cp -f ${CODEPATH}/*.o ${WDPATH}/' >> ${pbs_file}
echo 'cp -f ${CODEPATH}/vesemd ${WDPATH}/' >> ${pbs_file}

echo >> ${pbs_file}

COUNTER=0
for we in 1.0 10.0 100.0
do

outpath=output/${output_dir}/${model}${prob_num}/we${we}

for N in 4 8 12 16
do

for K in K2_Kx1Ky2
do

for L in _L8 _L16 _L32 _L64
do

infile=${inpath}/${K}${L}.dat

for t in 0.1 0.01 0.001 0.0001
do

COUNTER=$(( COUNTER + 1 ))

echo outpath${COUNTER}='${CODEPATH}'/${outpath} >> ${pbs_file}
echo outfile${COUNTER}='${CODEPATH}'/${outpath}/${model}${prob_num}_we${we}_dt${t}_${K}${L}_N${N} >> ${pbs_file}
echo input_arg${COUNTER}=\'-input ${infile} -output output${COUNTER} -prob ${prob_num} -N ${N} -alphaz 1.0 -deltat ${t} -re ${re} -we ${we} -beta ${beta}\' >> ${pbs_file}

done
done
done
done
done

echo >> ${pbs_file}

# Change to scratch directory
echo 'cd ${WDPATH}' >> ${pbs_file}

echo >> ${pbs_file}

echo 'arg=$(eval echo \${input_arg${INDEX}})' >> ${pbs_file}
echo './vesemd ${arg} > screenoutput${INDEX}.out' >> ${pbs_file}

echo >> ${pbs_file}

echo 'echo Finished ${INDEX} with argument ${arg}' >> ${pbs_file}

echo >> ${pbs_file}
# Copy output into home directory
echo 'folder_check=$(eval echo \${outpath${INDEX}})' >> ${pbs_file}
echo 'if [ ! -d ${folder_check} ]; then' >> ${pbs_file}
echo 'mkdir -p ${folder_check}' >> ${pbs_file}
echo 'fi' >> ${pbs_file}

echo >> ${pbs_file}

echo 'out=$(eval echo \${outfile${INDEX}})' >> ${pbs_file}
echo 'cp output${INDEX}_tec.dat ${out}_tec.dat' >> ${pbs_file}
echo 'cp output${INDEX}_tec_fine.dat ${out}_tec_fine.dat' >> ${pbs_file}
echo 'cp screenoutput${INDEX}.out ${out}.out' >> ${pbs_file}

echo >> ${pbs_file}
# Remove working directory in scratch.
echo 'rm -rf $WDPATH' >> ${pbs_file}


