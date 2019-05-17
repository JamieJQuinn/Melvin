#!/usr/bin/env bash

save_folder="data/salt_fingering_unstable"

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.js
{
  "Pr":1,
  "Ra":1.1e2,
  "RaXi":1e3,
  "tau":1e-1,
  "aspectRatio":1.41421356237,
  "initialDt":3e-6,

  "nN":51,
  "nZ":101,

  "icFile":"$save_folder/ICn1nZ101nN51_SF",
  "saveFolder":"$save_folder/",

  "timeBetweenSaves":0.1,
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
{ time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log

echo "==================== Comparing results"
cat << EOF > $save_folder/benchmark.txt
0   3.20258E-01   0.00000E+00   0.00000E+00
1   1.67806E-02   -7.18384E+00  -3.63938E-01
2   -1.11428E-07  1.97815E-04   2.17904E-06
3   1.10522E-09   -3.02042E-06  -3.61897E-08
4   7.32543E-12   -6.23672E-08  -2.83749E-10
5   -8.98958E-14  1.02699E-09   4.19486E-12
6   -7.70814E-16  7.38293E-12   1.68919E-14
7   -2.46911E-17  3.34380E-13   6.76961E-16
8   -1.31778E-18  3.01811E-14   4.21929E-17
9   -6.08775E-20  1.13005E-15   1.52551E-18
10  -2.25129E-21  6.70094E-17   6.59763E-20
11  3.43816E-21   -2.60157E-17  -4.59234E-20
12  8.36862E-20   4.02492E-16   -8.71710E-19
13  1.80054E-18   1.04760E-14   -1.81177E-17
14  3.59578E-17   -4.78009E-13  -4.81594E-16
15  6.26592E-16   -2.12341E-11  -1.01536E-14
16  9.13372E-15   -2.65792E-10  -1.28517E-13
17  1.19745E-13   1.66272E-09   -8.08902E-13
18  1.72416E-12   7.16886E-08   -8.86847E-12
19  2.82313E-11   -1.92940E-07  -4.49379E-10
20  4.19234E-10   -3.51243E-05  -1.19393E-08
EOF

comparison_results=$(python3 tools/print_variables.py $save_folder/dump0005.dat --max_print_mode 20 --n_modes 51 --n_gridpoints 101 | column -t | diff - $save_folder/benchmark.txt)

if [ "$comparison_results" ]; then
  echo "Simulation returned different results!"
  exit -1
else
  echo "Results match!"
  exit 0
fi
