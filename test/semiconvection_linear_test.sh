#!/usr/bin/env bash

set -e

save_folder="test/benchmark"

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.js
{
  "Pr":0.5,
  "Ra":3000,
  "aspectRatio":3,
  "icFile":"initial_conditions/ICn1nZ101nN51",
  "initialDt":1e-5,
  "nN":51,
  "nZ":101,
  "saveFolder":"test/benchmark/",
  "timeBetweenSaves":0.01,
  "totalTime":10,
  "xiGrad":-1,
  "tempGrad":-1,
  "RaXi":3e5,
  "tau":0.01
}
EOF

constants_file=$save_folder/constants.js

echo "==================== Building program"
make ddcLinear

echo "==================== Starting program"
time build/exe --constants $constants_file | tee $save_folder/log
