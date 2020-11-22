#!/usr/bin/env bash

save_folder="data/double_periodic_test"
constants_file=$save_folder/constants.json
nN=85
nZ=256
nX=256

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.json
{
  "Pr":7,
  "Ra":1.3e5,
  "RaXi":1e7,
  "tau":1e-2,
  "aspectRatio":0.5,
  "initialDt":1e-7,

  "nN":$nN,
  "nX":$nX,
  "nZ":$nZ,

  "verticalBoundaryConditions":"periodic",
  "horizontalBoundaryConditions":"periodic",
  "temperatureGradient":1,
  "salinityGradient":1,

  "icFile":"$save_folder/initial_conditions.dat",
  "saveFolder":"$save_folder/",

  "timeBetweenSaves":0.001,
  "totalTime":0.1,

  "isNonlinear":true,
  "isDoubleDiffusion":true
}
EOF

python tools/make_initial_conditions.py --output $save_folder/initial_conditions.dat --n_modes $nN --n_gridpoints $nZ --modes 50 --amp 0.01 --salt_fingering --periodic

echo "==================== Building program"
make -j4 release

echo "==================== Starting program"
{ time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log.txt
