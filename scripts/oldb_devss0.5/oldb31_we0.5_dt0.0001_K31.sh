# Script for prob31 we=0.5 dt=0.0001 K=31
export OMP_NUM_THREADS=1
cd ../../

./vesemd -input input/fixedcylinder/fixedcylindermeshK31.dat -output output/oldb_devss0.5/oldb31/we0.5/oldb31_we0.5_dt0.0001_K31_N13 -prob 31 -N 13 -alphaz 1.0 -deltat 0.0001 -re 0.0 -we 0.5 -beta 0.59 -beta_s 0.5 > output/oldb_devss0.5/oldb31/we0.5/oldb31_we0.5_dt0.0001_K31_N13.out && echo Finished We=0.5 N=13 K=31 Deltat=0.0001 &
./vesemd -input input/fixedcylinder/fixedcylindermeshK31.dat -output output/oldb_devss0.5/oldb31/we0.5/oldb31_we0.5_dt0.0001_K31_N14 -prob 31 -N 14 -alphaz 1.0 -deltat 0.0001 -re 0.0 -we 0.5 -beta 0.59 -beta_s 0.5 > output/oldb_devss0.5/oldb31/we0.5/oldb31_we0.5_dt0.0001_K31_N14.out && echo Finished We=0.5 N=14 K=31 Deltat=0.0001 &
./vesemd -input input/fixedcylinder/fixedcylindermeshK31.dat -output output/oldb_devss0.5/oldb31/we0.5/oldb31_we0.5_dt0.0001_K31_N15 -prob 31 -N 15 -alphaz 1.0 -deltat 0.0001 -re 0.0 -we 0.5 -beta 0.59 -beta_s 0.5 > output/oldb_devss0.5/oldb31/we0.5/oldb31_we0.5_dt0.0001_K31_N15.out && echo Finished We=0.5 N=15 K=31 Deltat=0.0001 &
./vesemd -input input/fixedcylinder/fixedcylindermeshK31.dat -output output/oldb_devss0.5/oldb31/we0.5/oldb31_we0.5_dt0.0001_K31_N16 -prob 31 -N 16 -alphaz 1.0 -deltat 0.0001 -re 0.0 -we 0.5 -beta 0.59 -beta_s 0.5 > output/oldb_devss0.5/oldb31/we0.5/oldb31_we0.5_dt0.0001_K31_N16.out && echo Finished We=0.5 N=16 K=31 Deltat=0.0001 &
wait
