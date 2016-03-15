#!/bin/bash
# Rayleigh numbers above 1e7 give horrible data
for i in `seq 2 6`;
do 
	folder=Ra1e$i
	mkdir -p $folder
	rm -f $folder/KineticEnergy.dat
	./box 101 51 3e-6 1e$i 0.5 3 3 $folder/ | tee $folder/log
done
# Need smalled time step for 
i=7
folder=Ra1e$i
mkdir -p $folder
rm -f $folder/KineticEnergy.dat
./box 101 51 5e-7 1e$i 0.5 3 3 $folder/ | tee $folder/log
