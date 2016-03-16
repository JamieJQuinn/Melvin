#/bin/bash
nZ=101
nN=51
dt=3e-6
Ra=1e6
Pr=0.5
a=3
T=3e-3 # Run for a thousand time steps
S=100 # Don't save
folder=TestData/
time ./box "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-o" $folder
