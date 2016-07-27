for file in ./*23*.pbs
do
qsub ${file}
done
