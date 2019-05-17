#!/usr/bin/env bash

save_folder="data/salt_fingering_unstable"

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.js
{
  "Pr":1,
  "Ra":1.1e4,
  "RaXi":1e6,
  "tau":1e-2,
  "aspectRatio":1.41421356237,
  "initialDt":3e-6,

  "nN":51,
  "nZ":101,

  "icFile":"$save_folder/ICn1nZ101nN51_SF",
  "saveFolder":"$save_folder/",

  "timeBetweenSaves":0.01,
  "totalTime":10,

  "isNonlinear":true,
  "isDoubleDiffusion":true
}
EOF

constants_file=$save_folder/constants.js
python tools/make_initial_conditions.py --output $save_folder/ICn1nZ101nN51_SF --n_modes 51 --n_gridpoints 101 --modes 1 25 --salt_fingering

echo "==================== Building program"
make clean
make release

echo "==================== Starting program"
{ time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log.txt
