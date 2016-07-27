for file in ./*21*.pbs
do
qsub ${file}
done
