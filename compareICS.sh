#/bin/bash
nZ=101
nN=51
dt=3e-6
Pr=0.5
a=3
T=1.5
S=0.05
Ra=1e6
for IC in ICn1nZ101nN51; do
		folder=n1Stable
		mkdir -p $folder
		rm -f $folder/*
		./box "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-i" $IC "-o" $folder/ | tee $folder/log
done
