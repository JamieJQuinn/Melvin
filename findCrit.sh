#/bin/bash
nZ=101
nN=51
dt=1e-5
Pr=0.5
a=3
T=10
S=0.1
for ICs in ICn1_2nZ101nN51 ICn1_8nZ101nN51 ; do
for i in 600 900 1300; do
	Ra="$i"
	folder=findCrit/"$ICs"Ra"$Ra"/
	mkdir -p $folder
	rm -f $folder/*
	./box "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-i" $ICs "-o" $folder/ | tee $folder/log
done
done
