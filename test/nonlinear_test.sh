#!/usr/bin/env bash

set -e

save_folder="test/benchmark"

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.js
{
  "Pr":0.5,
  "Ra":1000000,
  "aspectRatio":3,
  "icFile":"initial_conditions/ICn1nZ101nN51",
  "initialDt":3e-06,
  "nN":51,
  "nZ":101,
  "saveFolder":"test/benchmark/",
  "timeBetweenSaves":0.01,
  "totalTime":0.05,
  "tempGrad":-1
}
EOF

constants_file=$save_folder/constants.js

echo "==================== Building program"
make clean
make nonlinear

echo "==================== Starting program"
time build/exe --constants $constants_file | tee $save_folder/log

echo "==================== Comparing results"
comparison_results=$(python tools/print_variables.py test/benchmark/vars5.dat --max_print_mode 20 --n_modes 51 --n_gridpoints 101 | column -t | diff -q - test/benchmark_t5.txt)

if [ -n "$comparison_results" ]; then
  echo "Simulation returned different results!"
  exit -1
else
  echo "Results match!"
  exit 0
fi
