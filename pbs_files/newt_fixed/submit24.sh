for file in ./*24*.pbs
do
qsub ${file}
done
