#!/usr/bin/env bash
nZ=101
nN=51
dt=3e-6
Pr=0.5
a=3
T=3
S=0.01
Ra=1e6
folder=benchmark
ICs=initial_conditions/ICn1nZ101nN51
mkdir -p $folder
rm -f $folder/*
echo "Starting program"
build/exe "-nZ" $nZ "-nN" $nN "-dt" $dt "-Ra" $Ra "-Pr" $Pr "-a" $a "-T" $T "-S" $S "-i" $ICs "-o" $folder/ | tee $folder/log
