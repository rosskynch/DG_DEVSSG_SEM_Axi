# Script for prob31 we=0.6 dt=0.0001 K=26
export OMP_NUM_THREADS=1
cd ../../

./vesemd -input input/fixedcylinder/fixedcylindermeshK26.dat -output output/oldb_devss0.5/oldb31/we0.6/oldb31_we0.6_dt0.0001_K26_N15 -prob 31 -N 15 -alphaz 1.0 -deltat 0.0001 -re 0.0 -we 0.6 -beta 0.59 -beta_s 0.5 > output/oldb_devss0.5/oldb31/we0.6/oldb31_we0.6_dt0.0001_K26_N15.out && echo Finished We=0.6 N=15 K=26 Deltat=0.0001 &
./vesemd -input input/fixedcylinder/fixedcylindermeshK26.dat -output output/oldb_devss0.5/oldb31/we0.6/oldb31_we0.6_dt0.0001_K26_N16 -prob 31 -N 16 -alphaz 1.0 -deltat 0.0001 -re 0.0 -we 0.6 -beta 0.59 -beta_s 0.5 > output/oldb_devss0.5/oldb31/we0.6/oldb31_we0.6_dt0.0001_K26_N16.out && echo Finished We=0.6 N=16 K=26 Deltat=0.0001 &
./vesemd -input input/fixedcylinder/fixedcylindermeshK26.dat -output output/oldb_devss0.5/oldb31/we0.6/oldb31_we0.6_dt0.0001_K26_N17 -prob 31 -N 17 -alphaz 1.0 -deltat 0.0001 -re 0.0 -we 0.6 -beta 0.59 -beta_s 0.5 > output/oldb_devss0.5/oldb31/we0.6/oldb31_we0.6_dt0.0001_K26_N17.out && echo Finished We=0.6 N=17 K=26 Deltat=0.0001 &
./vesemd -input input/fixedcylinder/fixedcylindermeshK26.dat -output output/oldb_devss0.5/oldb31/we0.6/oldb31_we0.6_dt0.0001_K26_N18 -prob 31 -N 18 -alphaz 1.0 -deltat 0.0001 -re 0.0 -we 0.6 -beta 0.59 -beta_s 0.5 > output/oldb_devss0.5/oldb31/we0.6/oldb31_we0.6_dt0.0001_K26_N18.out && echo Finished We=0.6 N=18 K=26 Deltat=0.0001 &
wait
