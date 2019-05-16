#!/usr/bin/env bash

set -e

save_folder="test/benchmark"

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.js
{
  "Pr":0.5,
  "Ra":1,
  "aspectRatio":3,
  "icFile":"$save_folder/ICn1nZ101nN51_SF",
  "initialDt":1e-5,
  "nN":51,
  "nZ":101,
  "saveFolder":"test/benchmark/",
  "timeBetweenSaves":0.01,
  "totalTime":10,
  "RaXi":1000,
  "tau":0.1
}
EOF

constants_file=$save_folder/constants.js
python tools/make_initial_conditions.py --output $save_folder/ICn1nZ101nN51_SF --n_modes 51 --n_gridpoints 101 --salt_fingering --linear_stability

echo "==================== Building program"
make clean
make ddcLinear

echo "==================== Starting program"
time build/exe --constants $constants_file | tee $save_folder/log
