#!/usr/bin/env bash

save_folder="test/benchmark"
constants_file=$save_folder/constants.json
nN=50
nZ=101

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $constants_file
{
  "Pr":1.0,
  "Ra":1e6,
  "aspectRatio":2.0,
  "icFile":"$save_folder/initial_conditions.dat",
  "initialDt":3e-6,
  "nN":$nN,
  "nZ":$nZ,
  "saveFolder":"$save_folder/",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":0.05,
  "horizontalBoundaryConditions":"periodic"
}
EOF

python3 tools/make_initial_conditions.py --output $save_folder/initial_conditions.dat --perturb_vorticity --n_modes $nN --n_gridpoints $nZ --modes 1

echo "==================== Building program"
#make clean
make -j4 release

echo "==================== Starting program"
{ /usr/bin/time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log

echo "==================== Comparing results"

cat << EOF > $save_folder/benchmark.txt
0   6.33942E-01   -2.58680E+03  -2.38708E+02
1   -4.88682E-02  -1.88236E+02  -1.51581E+00
2   1.56070E-02   -6.31091E+01  -2.55162E+00
3   2.04740E-03   1.24106E+02   -4.76251E-03
4   -6.50115E-03  -6.28768E+01  1.66072E-01
5   1.99261E-03   -3.15210E+01  -1.08137E-02
6   1.63381E-03   4.81532E+01   -9.33626E-03
7   -1.36141E-03  -4.60071E+00  5.47044E-03
8   -1.76644E-04  -2.09012E+01  -2.31994E-03
9   5.88428E-04   8.90307E+00   -1.03457E-03
10  -9.22712E-05  6.57862E+00   1.49148E-03
11  -2.07816E-04  -5.54061E+00  -6.98422E-05
12  8.34180E-05   -1.43729E+00  -5.43916E-04
13  6.30868E-05   2.59056E+00   1.54545E-04
14  -4.46724E-05  6.74328E-02   1.58557E-04
15  -1.58322E-05  -1.05161E+00  -8.50090E-05
16  1.98845E-05   1.56931E-01   -3.82378E-05
17  2.52059E-06   3.88604E-01   3.54930E-05
18  -7.90862E-06  -1.25050E-01  6.73658E-06
19  3.99081E-07   -1.30612E-01  -1.29241E-05
20  2.85946E-06   6.91391E-02   -1.04342E-07
EOF

comparison_results=$(python3 tools/print_variables.py $save_folder/dump0005.dat --max_print_mode 20 --constants $save_folder/constants.json | column -t | diff - $save_folder/benchmark.txt)

if [ "$comparison_results" ]; then
  echo "Simulation returned different results!"
  exit -1
else
  echo "Results match!"
  exit 0
fi
