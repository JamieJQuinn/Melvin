#!/bin/bash
for i in `seq 4 10`;
do 
	folder=Ra1e$i
	mkdir -p $folder
	./box 101 51 3e-6 1e$i 0.5 3 3 $folder/ | tee $folder/log
done
