prob_num=24
for re in 0.001 0.01 0.1 1.0 10.0
do
echo re=${re} > newt${prob_num}_re${re}.sh
echo prob_num=${prob_num} >> newt${prob_num}_re${re}.sh
echo inpath=input/3dpipe/ >> newt${prob_num}_re${re}.sh
echo outpath=output/newt_fixed/newt${prob_num}/re${re}/ >> newt${prob_num}_re${re}.sh
cat prepfile23_24.txt >> newt${prob_num}_re${re}.sh
done


