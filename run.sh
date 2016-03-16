#/bin/bash
nZ=101
nN=51
dt=1e-7
Ra=1e7
Pr=0.5
a=3
T=3
S=0.1
folder=n1Turb/
ICs=n1StablenN101nZ101
mkdir -p $folder
rm -f $folder/*
echo "Starting program"
./box "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-i" $ICs "-o" $folder/ | tee $folder/log
