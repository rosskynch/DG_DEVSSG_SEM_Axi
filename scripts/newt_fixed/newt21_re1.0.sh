re=1.0
prob_num=21
inpath=input/fixedcylinder/fixedcylindermeshK
outpath=output/newt_fixed/newt21/re1.0/

export OMP_NUM_THREADS=1

# Now use alphaz=1.0
# move to parent dir.
cd ../../

for N in 4 8 12 16
do
for K in 7 10 20 26
do
for t in 0.1 0.01 0.001 0.0001
do

outfile_screen=${outpath}newt${prob_num}_Re${re}_dt${t}_K${K}_N${N}.out
outfile=${outpath}newt${prob_num}_Re${re}_dt${t}_K${K}_N${N}
infile=${inpath}${K}.dat

./vesemd -input ${infile} -output ${outfile} -prob ${prob_num} -N ${N} -alphaz 1.0 -deltat ${t} -re ${re} > ${outfile_screen}
done
done
done
