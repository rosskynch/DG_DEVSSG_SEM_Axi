mkdir ../output/newt_fixed/
for prob_num in 21 22 23 24
do
mkdir ../output/newt_fixed/newt${prob_num}/
for re in 0.001 0.01 0.1 1.0 10.0
do
mkdir ../output/newt_fixed/newt${prob_num}/re${re}/
done
done

