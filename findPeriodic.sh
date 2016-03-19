#/bin/bash
nZ=101
nN=101
dt=5e-8
Pr=0.5
a=3
T=0.8
S=0.01
for i in 7.5 7.75 8; do
	Ra="$i"e7
	folder=n1TurbRa"$Ra"/
	ICs=n1StablenN101nZ101 
	mkdir -p $folder
	rm -f $folder/*
	echo "Starting program"
	./box "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-i" $ICs "-o" $folder/ | tee $folder/log
done
