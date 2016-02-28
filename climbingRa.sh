#!/bin/bash
for i in `seq 0 10`;
do 
	./box 101 51 3e-6 1e$i 0.5 3 3 > $i.log
done
