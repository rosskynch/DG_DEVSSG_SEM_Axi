for file in ./*22*.pbs
do
qsub ${file}
done
