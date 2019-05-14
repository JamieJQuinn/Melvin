#!/usr/bin/env bash

set -e

nZ=101
nN=51
dt=3e-6
Pr=0.5
a=3
T=5e-2
S=0.01
Ra=1e6
folder=test/benchmark
ICs=initial_conditions/ICn1nZ101nN51
mkdir -p $folder
rm -f $folder/*
echo "==================== Building program"
make nonlinear
echo "==================== Starting program"
time build/exe -nZ $nZ -nN $nN -dt $dt -Ra $Ra -Pr $Pr -a $a -T $T -S $S -i $ICs -o $folder/ | tee $folder/log
echo "==================== Comparing results"
comparison_results=$(python tools/print_variables.py test/benchmark/vars5.dat --max_print_mode 20 --n_modes 51 --n_gridpoints 101 | column -t | diff -q - test/benchmark_t5.txt)
if [ -n "$comparison_results" ]; then
  echo "Simulation returned different results!"
  exit -1
else
  echo "Results match!"
  exit 0
fi
