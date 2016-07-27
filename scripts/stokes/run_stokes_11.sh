export OMP_NUM_THREADS=1
#move to parent dir.
cd ../../

for a in 0.0 0.001 0.01 0.1 1.0
do
#make sure output file is fresh.
rm output/stokes_11/alphaz_${a}/stokes_11_alphaz_${a}.out

for K in 7 10 20 26 31
do
for N in 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
./vesemd -input input/fixedcylinder/fixedcylindermeshK${K}.dat -prob 11 -output output/stokes_11/alphaz_${a}/stokes_11_K${K}_N${N} -N ${N} -alphaz ${a} >> output/stokes_11/alphaz_${a}/stokes_11_alphaz_${a}.out
done
done
done
