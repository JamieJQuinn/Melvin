#/bin/bash
nZ=101
nN=51
dt=1e-5
Pr=0.5
Ra=40000
RaXi=2e4
tau=1e-4
a=3
T=100
S=100
ICs=ICn1
folder=CriticalDCC/
mkdir -p $folder
rm -f $folder/*
./box_ddc "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-RaXi" $RaXi "-tau" $tau "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-i" $ICs "-o" $folder/ | tee $folder/log
