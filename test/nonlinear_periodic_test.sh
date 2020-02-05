#!/usr/bin/env bash

save_folder="test/benchmark"
nN=51
nZ=101

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.json
{
  "Pr":0.5,
  "Ra":1000000,
  "aspectRatio":3,
  "icFile":"$save_folder/ICn1nZ101nN51",
  "initialDt":3e-06,
  "nN":$nN,
  "nZ":$nZ,
  "saveFolder":"$save_folder/",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":0.05,
  "boundaryConditions":"periodic",
  "temperatureGradient":-1
}
EOF

constants_file=$save_folder/constants.json
python3 tools/make_initial_conditions.py --output $save_folder/ICn1nZ101nN51 --n_modes $nN --n_gridpoints $nZ --modes 2 --periodic

echo "==================== Building program"
#make clean
make release

echo "==================== Starting program"
{ /usr/bin/time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log

echo "==================== Comparing results"

cat << EOF > $save_folder/benchmark.txt
0   0.00000E+00   0.00000E+00   0.00000E+00
1   0.00000E+00   0.00000E+00   0.00000E+00
2   5.14185E-03   2.83664E+02   6.52475E+01
3   0.00000E+00   0.00000E+00   0.00000E+00
4   5.30221E-06   7.64586E-01   1.41952E-02
5   0.00000E+00   0.00000E+00   0.00000E+00
6   4.94995E-06   8.50346E-01   1.04877E-02
7   0.00000E+00   0.00000E+00   0.00000E+00
8   2.56679E-06   4.22616E-01   3.88155E-03
9   0.00000E+00   0.00000E+00   0.00000E+00
10  2.95157E-06   7.92550E-01   3.92012E-03
11  0.00000E+00   0.00000E+00   0.00000E+00
12  -2.19193E-06  -1.40428E+00  -7.16413E-03
13  0.00000E+00   0.00000E+00   0.00000E+00
14  5.86574E-07   5.90185E-01   2.27012E-03
15  0.00000E+00   0.00000E+00   0.00000E+00
16  -1.09191E-06  -1.96699E+00  -6.16608E-03
17  0.00000E+00   0.00000E+00   0.00000E+00
18  -8.03028E-09  -3.37906E-01  -8.80628E-04
19  0.00000E+00   0.00000E+00   0.00000E+00
20  -4.43677E-08  -7.10153E-01  -1.50173E-03
EOF

comparison_results=$(python3 tools/print_variables.py $save_folder/dump0005.dat --max_print_mode 20 --constants $save_folder/constants.json | column -t | diff - $save_folder/benchmark.txt)

if [ "$comparison_results" ]; then
  echo "Simulation returned different results!"
  exit -1
else
  echo "Results match!"
  exit 0
fi
