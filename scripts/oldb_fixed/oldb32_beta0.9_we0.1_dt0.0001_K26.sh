# Script for prob32 we=0.1 dt=0.0001 K=26 N=16
export OMP_NUM_THREADS=1
cd ../../

./vesemd -input input/fixedsphere/fixedspheremeshK26.dat -output output/oldb_fixed/oldb32/beta0.9/we0.1/oldb32_we0.1_dt0.0001_K26_N15 -prob 32 -N 15 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 0.1 -beta 0.888888888888888888888889 > output/oldb_fixed/oldb32/beta0.9/we0.1/oldb32_we0.1_dt0.0001_K26_N15.out && echo Finished  Beta=0.888888888888888888888889 We=0.1 N=15 K=26 Deltat=0.0001 &
./vesemd -input input/fixedsphere/fixedspheremeshK26.dat -output output/oldb_fixed/oldb32/beta0.9/we0.1/oldb32_we0.1_dt0.0001_K26_N16 -prob 32 -N 16 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 0.1 -beta 0.888888888888888888888889 > output/oldb_fixed/oldb32/beta0.9/we0.1/oldb32_we0.1_dt0.0001_K26_N16.out && echo Finished  Beta=0.888888888888888888888889 We=0.1 N=16 K=26 Deltat=0.0001 &
./vesemd -input input/fixedsphere/fixedspheremeshK26.dat -output output/oldb_fixed/oldb32/beta0.9/we0.1/oldb32_we0.1_dt0.0001_K26_N17 -prob 32 -N 17 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 0.1 -beta 0.888888888888888888888889 > output/oldb_fixed/oldb32/beta0.9/we0.1/oldb32_we0.1_dt0.0001_K26_N17.out && echo Finished  Beta=0.888888888888888888888889 We=0.1 N=17 K=26 Deltat=0.0001 &
./vesemd -input input/fixedsphere/fixedspheremeshK26.dat -output output/oldb_fixed/oldb32/beta0.9/we0.1/oldb32_we0.1_dt0.0001_K26_N18 -prob 32 -N 18 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 0.1 -beta 0.888888888888888888888889 > output/oldb_fixed/oldb32/beta0.9/we0.1/oldb32_we0.1_dt0.0001_K26_N18.out && echo Finished  Beta=0.888888888888888888888889 We=0.1 N=18 K=26 Deltat=0.0001 &
wait
