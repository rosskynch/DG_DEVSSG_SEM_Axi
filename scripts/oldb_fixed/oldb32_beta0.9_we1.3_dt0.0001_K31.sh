# Script for prob32 we=1.3 dt=0.0001 K=31 N=18
export OMP_NUM_THREADS=1
cd ../../

./vesemd -input input/fixedsphere/fixedspheremeshK31.dat -output output/oldb_fixed/oldb32/beta0.9/we1.3/oldb32_we1.3_dt0.0001_K31_N13 -prob 32 -N 13 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 1.3 -beta 0.888888888888888888888889 > output/oldb_fixed/oldb32/beta0.9/we1.3/oldb32_we1.3_dt0.0001_K31_N13.out && echo Finished  Beta=0.888888888888888888888889 We=1.3 N=13 K=31 Deltat=0.0001 &
./vesemd -input input/fixedsphere/fixedspheremeshK31.dat -output output/oldb_fixed/oldb32/beta0.9/we1.3/oldb32_we1.3_dt0.0001_K31_N14 -prob 32 -N 14 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 1.3 -beta 0.888888888888888888888889 > output/oldb_fixed/oldb32/beta0.9/we1.3/oldb32_we1.3_dt0.0001_K31_N14.out && echo Finished  Beta=0.888888888888888888888889 We=1.3 N=14 K=31 Deltat=0.0001 &
./vesemd -input input/fixedsphere/fixedspheremeshK31.dat -output output/oldb_fixed/oldb32/beta0.9/we1.3/oldb32_we1.3_dt0.0001_K31_N15 -prob 32 -N 15 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 1.3 -beta 0.888888888888888888888889 > output/oldb_fixed/oldb32/beta0.9/we1.3/oldb32_we1.3_dt0.0001_K31_N15.out && echo Finished  Beta=0.888888888888888888888889 We=1.3 N=15 K=31 Deltat=0.0001 &
./vesemd -input input/fixedsphere/fixedspheremeshK31.dat -output output/oldb_fixed/oldb32/beta0.9/we1.3/oldb32_we1.3_dt0.0001_K31_N16 -prob 32 -N 16 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 1.3 -beta 0.888888888888888888888889 > output/oldb_fixed/oldb32/beta0.9/we1.3/oldb32_we1.3_dt0.0001_K31_N16.out && echo Finished  Beta=0.888888888888888888888889 We=1.3 N=16 K=31 Deltat=0.0001 &
wait
