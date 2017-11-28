#/bin/bash
nZ=101
nN=101
dt=3e-6
Ra=1e6
Pr=0.5
a=3
T=3e-3 # Run for a thousand time steps
S=100 # Don't save
folder=TestData/
ICs=ICn1nZ101nN101
mkdir -p $folder
time ./box "-i" $ICs "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-o" $folder
