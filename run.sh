#/bin/bash
nZ=101
nN=51
dt=3e-6
Ra=1e5
Pr=0.5
a=3
T=0.4
S=0.002
folder=n1_8testCont/
ICs=n1_8InitData/20.dat
mkdir -p $folder
rm -f $folder/*
echo "Starting program"
./box "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-i" $ICs "-o" $folder/ | tee $folder/log
