#/bin/bash
nZ=101
nN=51
dt=3e-6
Ra=1e5
Pr=0.5
a=3
T=0.5
S=0.01
folder=n1/
ICs=ICn1
mkdir -p $folder
rm -f $folder/*
echo "Starting program"
./box "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-i" $ICs "-o" $folder/ | tee $folder/log
