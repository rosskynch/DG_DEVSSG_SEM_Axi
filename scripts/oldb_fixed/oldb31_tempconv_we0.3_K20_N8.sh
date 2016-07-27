#Temporal convergence script for prob31
export OMP_NUM_THREADS=1
cd ../../

./vesemd -input input/fixedcylinder/fixedcylindermeshK20.dat -output output/oldb_fixed/oldb31/tempconv/oldb31_we0.3_dt0.01_K20_N8 -prob 31 -N 8 -alphaz 1.0 -deltat 0.01 -re 0.0 -we 0.3 -beta 0.59 > output/oldb_fixed/oldb31/tempconv/oldb31_we0.3_dt0.01_K20_N8.out && echo Finished We=0.3 N=8 K=20 Deltat=0.01 &
./vesemd -input input/fixedcylinder/fixedcylindermeshK20.dat -output output/oldb_fixed/oldb31/tempconv/oldb31_we0.3_dt0.001_K20_N8 -prob 31 -N 8 -alphaz 1.0 -deltat 0.001 -re 0.0 -we 0.3 -beta 0.59 > output/oldb_fixed/oldb31/tempconv/oldb31_we0.3_dt0.001_K20_N8.out && echo Finished We=0.3 N=8 K=20 Deltat=0.001 &
./vesemd -input input/fixedcylinder/fixedcylindermeshK20.dat -output output/oldb_fixed/oldb31/tempconv/oldb31_we0.3_dt0.0001_K20_N8 -prob 31 -N 8 -alphaz 1.0 -deltat 0.0001 -re 0.0 -we 0.3 -beta 0.59 > output/oldb_fixed/oldb31/tempconv/oldb31_we0.3_dt0.0001_K20_N8.out && echo Finished We=0.3 N=8 K=20 Deltat=0.0001 &
./vesemd -input input/fixedcylinder/fixedcylindermeshK20.dat -output output/oldb_fixed/oldb31/tempconv/oldb31_we0.3_dt0.00001_K20_N8 -prob 31 -N 8 -alphaz 1.0 -deltat 0.00001 -re 0.0 -we 0.3 -beta 0.59 > output/oldb_fixed/oldb31/tempconv/oldb31_we0.3_dt0.00001_K20_N8.out && echo Finished We=0.3 N=8 K=20 Deltat=0.00001 &
wait
