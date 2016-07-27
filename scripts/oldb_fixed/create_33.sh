prob_num=33
CODEDIR=/home/sacrmk/FINAL
for N in 4 8 12 16
do
for L in _L8 _L16 _L32 _L64
do
for we in 1.0 10.0 100.0
do
echo we=${we} > oldb${prob_num}_we${we}${L}_N${N}.sh
echo prob_num=${prob_num} >> oldb${prob_num}_we${we}${L}_N${N}.sh
echo L=${L} >> oldb${prob_num}_we${we}${L}_N${N}.sh
echo N=${N} >> oldb${prob_num}_we${we}${L}_N${N}.sh
echo inpath=${CODEPATH}/input/2dchannel/ >> oldb${prob_num}_we${we}${L}_N${N}.sh
echo outpath=${CODEPATH}/output/oldb_fixed/oldb${prob_num}/we${we}/ >> oldb${prob_num}_we${we}${L}_N${N}.sh
cat prepfile33_34.txt >> oldb${prob_num}_we${we}${L}_N${N}.sh
done
done
done
