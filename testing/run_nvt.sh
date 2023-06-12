#!/bin/bash


T=400


for P in 2300; do
for n3 in 10 40 80 200 400 600 800; do

name=p-$p-n1-4000-n3-$N3-T-$T

natoms=`awk -v n3=$n3 'BEGIN{print n3 + 4000}'`
sed "s/Tx/$T/g; s/PX/$P/g; s/NatomsX/$natoms/g; s/StateNameX/$name/g; s/n3x/$n3/g" inputFile.template > inputFile-$name.json

./gr_GCMC_clion -i inputFile.test.json $n $box_size $T > log-$name.log &

done; done
wait
exit
