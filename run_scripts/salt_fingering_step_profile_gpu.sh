#!/usr/bin/env bash

save_folder="data/salt_fingering_step_profile_gpu"

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.js
{
  "Pr":1,
  "Ra":1e6,
  "RaXi":1e7,
  "tau":1e-2,
  "aspectRatio":1.41421356237,
  "initialDt":3e-6,

  "nN":128,
  "nZ":256,

  "icFile":"$save_folder/ICn1nZ101nN51_SF",
  "saveFolder":"$save_folder/",

  "timeBetweenSaves":0.001,
  "totalTime":1,

  "isNonlinear":true,
  "isDoubleDiffusion":true,

  "isCudaEnabled":true,
  "threadsPerBlock_x":16,
  "threadsPerBlock_y":32
}
EOF

constants_file=$save_folder/constants.js
python tools/make_initial_conditions.py --output $save_folder/ICn1nZ101nN51_SF --n_modes 101 --n_gridpoints 201 --modes 20 25 --amp 0.01 --salt_fingering --step_profile

echo "==================== Building program"
make gpu

echo "==================== Starting program"
{ time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log.txt
