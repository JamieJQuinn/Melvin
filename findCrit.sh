#/bin/bash
nZ=101
nN=101
dt=3e-6
Pr=0.5
a=3
T=1
S=0.01
for ICs in ICn1_8 ; do
for i in 550 750 1150 1250; do
	Ra="$i"
	folder=findCrit/"$ICs"Ra"$Ra"/
	mkdir -p $folder
	rm -f $folder/*
	./box "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-i" $ICs "-o" $folder/ | tee $folder/log
done
done
