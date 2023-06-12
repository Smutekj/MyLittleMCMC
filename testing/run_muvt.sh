#!/bin/bash


T=400

for L in 40  ; do
for mu in 0.01 2.0 ; do

name=L-$L-mu-$mu-L-$L
natoms=2000
sed "s/Tx/$T/g; s/Lx/$L/g; s/mux/$mu/g; s/natomsx/$natoms/g" inputFile.template > inputFile-$name.json

./gr_GCMC_clion -i inputFile-$name.json > log-$name.log &

done; done
wait
exit
