#/bin/bash
nZ=101
nN=51
dt=1e-5
Pr=0.5
a=3
T=100
S=100
Ra=16000
folder=criticalThermal
ICs=ICn1
mkdir -p $folder
rm -f $folder/*
./box_linear "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-i" $ICs "-o" $folder/ | tee $folder/log
