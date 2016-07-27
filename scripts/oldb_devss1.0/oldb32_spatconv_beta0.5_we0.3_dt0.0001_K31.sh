#Spatial convergence script for prob32
export OMP_NUM_THREADS=1
cd ../../

./vesemd -input input/fixedsphere/fixedspheremeshK31.dat -output output/oldb_devss1.0/oldb32/beta0.5/spatconv/oldb32_we0.3_dt0.0001_K31_N4 -prob 32 -N 4 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 0.3 -beta 0.5 -beta_s 1.0 > output/oldb_devss1.0/oldb32/beta0.5/spatconv/oldb32_we0.3_dt0.0001_K31_N4.out && echo Finished  Beta=0.5 We=0.3 N=4 K=31 Deltat=0.0001 &
./vesemd -input input/fixedsphere/fixedspheremeshK31.dat -output output/oldb_devss1.0/oldb32/beta0.5/spatconv/oldb32_we0.3_dt0.0001_K31_N8 -prob 32 -N 8 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 0.3 -beta 0.5 -beta_s 1.0 > output/oldb_devss1.0/oldb32/beta0.5/spatconv/oldb32_we0.3_dt0.0001_K31_N8.out && echo Finished  Beta=0.5 We=0.3 N=8 K=31 Deltat=0.0001 &
./vesemd -input input/fixedsphere/fixedspheremeshK31.dat -output output/oldb_devss1.0/oldb32/beta0.5/spatconv/oldb32_we0.3_dt0.0001_K31_N12 -prob 32 -N 12 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 0.3 -beta 0.5 -beta_s 1.0 > output/oldb_devss1.0/oldb32/beta0.5/spatconv/oldb32_we0.3_dt0.0001_K31_N12.out && echo Finished  Beta=0.5 We=0.3 N=12 K=31 Deltat=0.0001 &
./vesemd -input input/fixedsphere/fixedspheremeshK31.dat -output output/oldb_devss1.0/oldb32/beta0.5/spatconv/oldb32_we0.3_dt0.0001_K31_N16 -prob 32 -N 16 -alphaz 1.0 -deltat 0.0001 -re 0.01 -we 0.3 -beta 0.5 -beta_s 1.0 > output/oldb_devss1.0/oldb32/beta0.5/spatconv/oldb32_we0.3_dt0.0001_K31_N16.out && echo Finished  Beta=0.5 We=0.3 N=16 K=31 Deltat=0.0001 &
wait
