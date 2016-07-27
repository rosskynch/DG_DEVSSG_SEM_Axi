we=1.0
prob_num=34
L=_L8
N=16
inpath=input/3dpipe/
outpath=output/oldb_devss0.5/oldb34/we1.0/

export OMP_NUM_THREADS=1

# Now use alphaz=1.0
# move to parent dir.
cd ../../

re=1.0
beta=0.1111111111111111111111111111111111111111
beta_s=0.5

for K in K1 K2_Kx1Ky2 K2_Kx2Ky1 K4_Kx2Ky2
do
for t in 0.1 0.01 0.001 0.0001
do

outfile_screen=${outpath}oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}.out
outfile=${outpath}oldb${prob_num}_we${we}_dt${t}_${K}${L}_N${N}
infile=${inpath}${K}${L}.dat

./vesemd -input ${infile} -output ${outfile} -prob ${prob_num} -N ${N} -alphaz 1.0 -deltat ${t} -re ${re} -we ${we} -beta ${beta} -beta_s ${beta_s} > ${outfile_screen}

echo Finished We=${we} N=${N} ${K}${L} Deltat=${t}
done
done
