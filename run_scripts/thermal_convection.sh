#!/usr/bin/env bash

save_folder="data/thermal_convection"

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.json
{
  "Pr":0.5,
  "Ra":1e6,
  "aspectRatio":3,
  "icFile":"$save_folder/initial_conditions.dat",
  "initialDt":3e-06,
  "nN":51,
  "nZ":101,
  "saveFolder":"$save_folder/",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":1
}
EOF

constants_file=$save_folder/constants.json
python3 tools/make_initial_conditions.py --output $save_folder/initial_conditions.dat --n_modes 51 --n_gridpoints 101 --modes 1

echo "==================== Building program"
make clean
make release

echo "==================== Starting program"
{ /usr/bin/time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log
