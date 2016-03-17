#/bin/bash
nZ=101
nN=101
dt=5e-8
Pr=0.5
a=3
T=0.3
S=0.001
for i in `seq 1 2 10`; do
	Ra="$i"e7
	folder=n1TurbRa"$Ra"/
	ICs=n1StablenN101nZ101
	mkdir -p $folder
	rm -f $folder/*
	echo "Starting program"
	./box "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-i" $ICs "-o" $folder/ | tee $folder/log
done
