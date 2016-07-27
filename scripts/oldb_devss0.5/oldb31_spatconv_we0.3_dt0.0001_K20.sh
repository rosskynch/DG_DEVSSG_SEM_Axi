#Spatial convergence script for prob31
export OMP_NUM_THREADS=1
cd ../../

./vesemd -input input/fixedcylinder/fixedcylindermeshK20.dat -output output/oldb_devss0.5/oldb31/spatconv/oldb31_we0.3_dt0.0001_K20_N4 -prob 31 -N 4 -alphaz 1.0 -deltat 0.0001 -re 0.0 -we 0.3 -beta 0.59 -beta_s 0.5 > output/oldb_devss0.5/oldb31/spatconv/oldb31_we0.3_dt0.0001_K20_N4.out && echo Finished We=0.3 N=4 K=20 Deltat=0.0001 &
./vesemd -input input/fixedcylinder/fixedcylindermeshK20.dat -output output/oldb_devss0.5/oldb31/spatconv/oldb31_we0.3_dt0.0001_K20_N8 -prob 31 -N 8 -alphaz 1.0 -deltat 0.0001 -re 0.0 -we 0.3 -beta 0.59 -beta_s 0.5 > output/oldb_devss0.5/oldb31/spatconv/oldb31_we0.3_dt0.0001_K20_N8.out && echo Finished We=0.3 N=8 K=20 Deltat=0.0001 &
./vesemd -input input/fixedcylinder/fixedcylindermeshK20.dat -output output/oldb_devss0.5/oldb31/spatconv/oldb31_we0.3_dt0.0001_K20_N12 -prob 31 -N 12 -alphaz 1.0 -deltat 0.0001 -re 0.0 -we 0.3 -beta 0.59 -beta_s 0.5 > output/oldb_devss0.5/oldb31/spatconv/oldb31_we0.3_dt0.0001_K20_N12.out && echo Finished We=0.3 N=12 K=20 Deltat=0.0001 &
./vesemd -input input/fixedcylinder/fixedcylindermeshK20.dat -output output/oldb_devss0.5/oldb31/spatconv/oldb31_we0.3_dt0.0001_K20_N16 -prob 31 -N 16 -alphaz 1.0 -deltat 0.0001 -re 0.0 -we 0.3 -beta 0.59 -beta_s 0.5 > output/oldb_devss0.5/oldb31/spatconv/oldb31_we0.3_dt0.0001_K20_N16.out && echo Finished We=0.3 N=16 K=20 Deltat=0.0001 &
wait
