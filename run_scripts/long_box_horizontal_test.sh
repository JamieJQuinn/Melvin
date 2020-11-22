#!/usr/bin/env bash

# This is a test of a high aspect ratio box

save_folder="data/long_box_horizontal_test"
constants_file=$save_folder/constants.json
nN=340
nZ=512

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $constants_file
{
  "Pr":1.0,
  "Ra":1e8,
  "aspectRatio":10.0,
  "icFile":"$save_folder/initial_conditions.dat",
  "initialDt":1e-8,
  "nN":$nN,
  "nZ":$nZ,
  "nX":1024,
  "saveFolder":"$save_folder/",
  "timeBetweenSaves":1e-6,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":0.01,
  "horizontalBoundaryConditions":"impermeable"
}
EOF

python3 tools/make_initial_conditions.py --output $save_folder/initial_conditions.dat --n_modes $nN --n_gridpoints $nZ --modes 10

echo "==================== Building program"
#make clean
make -j4 release

echo "==================== Starting program"
{ /usr/bin/time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log
