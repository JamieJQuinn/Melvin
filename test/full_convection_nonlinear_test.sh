#!/usr/bin/env bash

set -e

save_folder="test/benchmark"

mkdir -p $save_folder
rm -f $save_folder/*

#cat << EOF > $save_folder/constants.js
#{
  #"Pr":1,
  #"Ra":5000,
  #"RaXi":10000,
  #"tau":0.1,
  #"aspectRatio":1.41421356237,
  #"icFile":"$save_folder/ICn1nZ101nN51_SF",
  #"initialDt":3e-06,
  #"nN":51,
  #"nZ":101,
  #"saveFolder":"test/benchmark/",
  #"timeBetweenSaves":0.1,
  #"totalTime":0.5
#}
#EOF

cat << EOF > $save_folder/constants.js
{
  "Pr":0.5,
  "Ra":50000,
  "RaXi":100000,
  "tau":0.01,
  "aspectRatio":3,
  "icFile":"$save_folder/ICn1nZ101nN51_SF",
  "initialDt":1e-5,
  "nN":51,
  "nZ":101,
  "saveFolder":"test/benchmark/",
  "timeBetweenSaves":0.1,
  "totalTime":10
}
EOF

constants_file=$save_folder/constants.js
python tools/make_initial_conditions.py --output $save_folder/ICn1nZ101nN51_SF --n_modes 51 --n_gridpoints 101 --modes 1 --combined_convection

echo "==================== Building program"
make clean
make ddcNonlinear

echo "==================== Starting program"
time build/exe --constants $constants_file | tee $save_folder/log

#echo "==================== Comparing results"
#comparison_results=$(python tools/print_variables.py test/benchmark/vars5.dat --max_print_mode 20 --n_modes 51 --n_gridpoints 101 | column -t | diff -q - test/benchmark_t5.txt)

#if [ -n "$comparison_results" ]; then
  #echo "Simulation returned different results!"
  #exit -1
#else
  #echo "Results match!"
  #exit 0
#fi
