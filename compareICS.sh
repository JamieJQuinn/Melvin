#/bin/bash
nZ=101
nN=51
dt=3e-6
Pr=0.5
a=3
T=10
S=5e-3
Ra=1e5
for IC in ICn1_8; do
	for nN in 51 101; do
		folder=ICsComparison/"$IC"nN"$nN"
		mkdir -p $folder
		rm -f $folder/*
		./box "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-i" $IC "-o" $folder/ | tee $folder/log
	done
done
