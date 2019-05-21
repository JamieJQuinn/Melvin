#!/usr/bin/env bash

save_folder="data/thermal_rayleigh_taylor"

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.js
{
  "Pr":0.5,
  "Ra":1e7,
  "aspectRatio":1.41,
  "initialDt":3e-6,

  "nN":51,
  "nZ":101,

  "icFile":"$save_folder/ICn1nZ101nN51",
  "saveFolder":"$save_folder/",

  "timeBetweenSaves":0.01,
  "totalTime":0.05,

  "isNonlinear":true,
  "isDoubleDiffusion":false
}
EOF

constants_file=$save_folder/constants.js
python tools/make_initial_conditions.py --output $save_folder/ICn1nZ101nN51 --n_modes 51 --n_gridpoints 101 --modes 10 15 --amp 0.1 --thermal_rayleigh_taylor

echo "==================== Building program"
#make clean
make release

echo "==================== Starting program"
{ time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log.txt
