#!/usr/bin/env bash

save_folder="test/benchmark"
constants_file=$save_folder/constants.json
nN=170
nZ=256

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $constants_file
{
  "Pr":1.0,
  "Ra":1e7,
  "aspectRatio":2.0,
  "icFile":"$save_folder/initial_conditions.dat",
  "initialDt":1e-7,
  "nN":$nN,
  "nZ":$nZ,
  "nX":512,
  "saveFolder":"$save_folder/",
  "timeBetweenSaves":0.0001,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":0.1,
  "horizontalBoundaryCondi--vorticitytions":"periodic"
}
EOF

python3 tools/make_initial_conditions.py --output $save_folder/initial_conditions.dat --perturb_vorticity --n_modes $nN --n_gridpoints $nZ --modes 1

echo "==================== Building program"
#make clean
make -j4 release

echo "==================== Starting program"
{ /usr/bin/time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log
