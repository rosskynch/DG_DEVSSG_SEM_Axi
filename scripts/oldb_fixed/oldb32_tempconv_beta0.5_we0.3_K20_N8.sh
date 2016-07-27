#Temporal convergence script for prob32
export OMP_NUM_THREADS=1
cd ../../

./vesemd -input input/fixedsphere/fixedspheremeshK20.dat -output output/oldb_fixed/oldb32/beta0.5/tempconv/oldb32_we0.3_dt0.01_K20_N8 -prob 32 -N 8 -alphaz 1.0 -deltat 0.01 -re 0.01 -we 0.3 -beta 0.5 > output/oldb_fixed/oldb32/beta0.5/tempconv/oldb32_we0.3_dt0.01_K20_N8.out && echo Finished  Beta=0.5 We=0.3 N=8 K=20 Deltat=0.01 &
./vesemd -input input/fixedsphere/fixedspheremeshK20.dat -output output/oldb_fixed/oldb32/beta0.5/tempconv/oldb32_we0.3_dt0.001_K20_N8 -prob 32 -N 8 -alphaz 1.0 -deltat 0.001 -re 0.01 -we 0.3 -beta 0.5 > output/oldb_fixed/oldb32/beta0.5/tempconv/oldb32_we0.3_dt0.001_K20_N8.out && echo Finished  Beta=0.5 We=0.3 N=8 K=20 Deltat=0.001 &
./vesemd -input input/fixedsphere/fixedspheremeshK20.dat -output output/oldb_fixed/oldb32/beta0.5/tempconv/oldb32_we0.3_dt0.0001_K20_N8 -prob 32 -N 8 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 0.3 -beta 0.5 > output/oldb_fixed/oldb32/beta0.5/tempconv/oldb32_we0.3_dt0.0001_K20_N8.out && echo Finished  Beta=0.5 We=0.3 N=8 K=20 Deltat=0.0001 &
./vesemd -input input/fixedsphere/fixedspheremeshK20.dat -output output/oldb_fixed/oldb32/beta0.5/tempconv/oldb32_we0.3_dt0.00001_K20_N8 -prob 32 -N 8 -alphaz 1.0 -deltat 0.00001 -re 0.01 -we 0.3 -beta 0.5 > output/oldb_fixed/oldb32/beta0.5/tempconv/oldb32_we0.3_dt0.00001_K20_N8.out && echo Finished  Beta=0.5 We=0.3 N=8 K=20 Deltat=0.00001 &
wait
