#!/bin/bash
for i in `seq 8 10`;
do 
	folder=Ra1e$i
	mkdir -p $folder
	./box 101 51 5e-7 1e$i 0.5 3 3 $folder/ | tee $folder/log
done
