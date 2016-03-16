#!/bin/bash
for i in `seq 600 100 700`;
do 
	Ra="$i"
	folder=n1_8/Ra$Ra
	mkdir -p $folder
	rm -f $folder/KineticEnergy.dat
	./box 101 51 3e-6 $Ra 0.5 3 2 $folder/ | tee $folder/log
done
# Need smalled time step for 
#i=7
#folder=Ra1e$i
#mkdir -p $folder
#rm -f $folder/KineticEnergy.dat
#./box 101 51 5e-7 1e$i 0.5 3 3 $folder/ | tee $folder/log
