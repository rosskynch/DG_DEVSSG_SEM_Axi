export OMP_NUM_THREADS=1
#move to parent dir.
cd ../../

for a in 0.0 0.001 0.01 0.1 1.0
do
#make sure output file is fresh.
rm output/stokes_14/alphaz_${a}/stokes_14_alphaz_${a}.out

for K in 5 8 10 16
do
for N in 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
./vesemd -input input/model_cyl/fixedcylindermeshK${K}.dat -prob 14 -output output/stokes_14/alphaz_${a}/stokes_14_K${K}_N${N} -N ${N} -alphaz ${a} >> output/stokes_14/alphaz_${a}/stokes_14_alphaz_${a}.out
done
done
done
