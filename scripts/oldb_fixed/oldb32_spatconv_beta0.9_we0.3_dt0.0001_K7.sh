#Spatial convergence script for prob32
export OMP_NUM_THREADS=1
cd ../../

./vesemd -input input/fixedsphere/fixedspheremeshK7.dat -output output/oldb_fixed/oldb32/beta0.9/spatconv/oldb32_we0.3_dt0.0001_K7_N4 -prob 32 -N 4 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 0.3 -beta 0.888888888888888888888889 > output/oldb_fixed/oldb32/beta0.9/spatconv/oldb32_we0.3_dt0.0001_K7_N4.out && echo Finished  Beta=0.888888888888888888888889 We=0.3 N=4 K=7 Deltat=0.0001 &
./vesemd -input input/fixedsphere/fixedspheremeshK7.dat -output output/oldb_fixed/oldb32/beta0.9/spatconv/oldb32_we0.3_dt0.0001_K7_N8 -prob 32 -N 8 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 0.3 -beta 0.888888888888888888888889 > output/oldb_fixed/oldb32/beta0.9/spatconv/oldb32_we0.3_dt0.0001_K7_N8.out && echo Finished  Beta=0.888888888888888888888889 We=0.3 N=8 K=7 Deltat=0.0001 &
./vesemd -input input/fixedsphere/fixedspheremeshK7.dat -output output/oldb_fixed/oldb32/beta0.9/spatconv/oldb32_we0.3_dt0.0001_K7_N12 -prob 32 -N 12 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 0.3 -beta 0.888888888888888888888889 > output/oldb_fixed/oldb32/beta0.9/spatconv/oldb32_we0.3_dt0.0001_K7_N12.out && echo Finished  Beta=0.888888888888888888888889 We=0.3 N=12 K=7 Deltat=0.0001 &
./vesemd -input input/fixedsphere/fixedspheremeshK7.dat -output output/oldb_fixed/oldb32/beta0.9/spatconv/oldb32_we0.3_dt0.0001_K7_N16 -prob 32 -N 16 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 0.3 -beta 0.888888888888888888888889 > output/oldb_fixed/oldb32/beta0.9/spatconv/oldb32_we0.3_dt0.0001_K7_N16.out && echo Finished  Beta=0.888888888888888888888889 We=0.3 N=16 K=7 Deltat=0.0001 &
wait
