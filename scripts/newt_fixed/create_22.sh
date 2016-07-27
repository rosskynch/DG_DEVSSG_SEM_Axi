prob_num=22
for re in 0.001 0.01 0.1 1.0 10.0
do
echo re=${re} > newt${prob_num}_re${re}.sh
echo prob_num=${prob_num} >> newt${prob_num}_re${re}.sh
echo inpath=input/fixedsphere/fixedspheremeshK >> newt${prob_num}_re${re}.sh
echo outpath=output/newt_fixed/newt${prob_num}/re${re}/ >> newt${prob_num}_re${re}.sh
cat prepfile21_22.txt >> newt${prob_num}_re${re}.sh
done

