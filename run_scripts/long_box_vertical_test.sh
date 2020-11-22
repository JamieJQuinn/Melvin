#!/usr/bin/env bash

# This is a test of a low aspect ratio box

save_folder="data/long_box_vertical_test"
constants_file=$save_folder/constants.json
nN=85
nZ=256

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $constants_file
{
  "Pr":1.0,
  "Ra":1e6,
  "aspectRatio":0.1,
  "icFile":"$save_folder/initial_conditions.dat",
  "initialDt":1e-7,
  "nN":$nN,
  "nZ":$nZ,
  "nX":256,
  "saveFolder":"$save_folder/",
  "timeBetweenSaves":1e-5,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":0.01,
  "horizontalBoundaryConditions":"periodic"
}
EOF

python3 tools/make_initial_conditions.py --output $save_folder/initial_conditions.dat --n_modes $nN --n_gridpoints $nZ --modes 1

echo "==================== Building program"
#make clean
make -j4 release

echo "==================== Starting program"
{ /usr/bin/time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log
