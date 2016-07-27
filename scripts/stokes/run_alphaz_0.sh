export OMP_NUM_THREADS=1
#move to parent dir.
cd ../../

for a in 0.0
do
rm output/stokes_11/alphaz_${a}/stokes_11_alphaz_${a}.out
rm output/stokes_12/alphaz_${a}/stokes_12_alphaz_${a}.out
rm output/stokes_13/alphaz_${a}/stokes_13_alphaz_${a}.out
rm output/stokes_14/alphaz_${a}/stokes_14_alphaz_${a}.out

for K in 7 10 20 26 31
do
for N in 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
./vesemd -input input/fixedcylinder/fixedcylindermeshK${K}.dat -prob 11 -output output/stokes_11/alphaz_${a}/stokes_11_K${K}_N${N} -N ${N} -alphaz ${a} >> output/stokes_11/alphaz_${a}/screen.out
done
done

for K in 7 10 20 26 31
do
for N in 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
./vesemd -input input/fixedsphere/fixedspheremeshK${K}.dat -prob 12 -output output/stokes_12/alphaz_${a}/stokes_12_K${K}_N${N} -N ${N} -alphaz ${a} >> output/stokes_12/alphaz_${a}/screen.out
done
done

for K in 1 2 4 16 16_nonuniform
do
for N in 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
./vesemd -input input/model_soln/unitsquareK${K}.dat -prob 13 -output output/stokes_13/alphaz_${a}/stokes_13_K${K}_N${N} -N ${N} -alphaz ${a} >> output/stokes_13/alphaz_${a}/screen.out
done
done

for K in 5 8 10 16
do
for N in 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
./vesemd -input input/model_cyl/fixedcylindermeshK${K}.dat -prob 14 -output output/stokes_14/alphaz_${a}/stokes_14_K${K}_N${N} -N ${N} -alphaz ${a} >> output/stokes_14/alphaz_${a}/screen.out
done
done

done
