#!/bin/bash
# Simulate for Ra = 2e6 -> 9e6
# Rayleigh numbers above 1e7 breach CFL, too slow
for j in `seq 6 7`;
do
for i in `seq 1 9`;
do 
	Ra="$i"e$j
	folder=squareRa$Ra
	mkdir -p $folder
	./box 101 51 3e-6 $Ra 1 1 1 $folder/ | tee $folder/log
done
done
#	folder=periodic
#	mkdir -p $folder
#	./box 101 51 1e-7 1e8 1 1 3 $folder/ | tee $folder/log
