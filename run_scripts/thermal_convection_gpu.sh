#!/usr/bin/env bash

save_folder="data/thermal_convection_gpu"

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.js
{
  "Pr":0.5,
  "Ra":1000000,
  "aspectRatio":3,
  "icFile":"$save_folder/ICn1nZ101nN51",
  "initialDt":3e-06,
  "nN":51,
  "nZ":101,
  "saveFolder":"$save_folder/",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":1,

  "isCudaEnabled":true,
  "threadsPerBlock_x":16,
  "threadsPerBlock_y":32
}
EOF

constants_file=$save_folder/constants.js
python3 tools/make_initial_conditions.py --output $save_folder/ICn1nZ101nN51 --n_modes 51 --n_gridpoints 101 --modes 1

echo "==================== Building program"
make clean
make gpu

echo "==================== Starting program"
{ /usr/bin/time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log
