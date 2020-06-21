#!/usr/bin/env bash

save_folder="test/benchmark"
constants_file=$save_folder/constants.json
ic_file=$save_folder/initial_conditions.dat
nN=50
nZ=101

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $constants_file
{
  "Pr":1.0,
  "Ra":1e6,
  "aspectRatio":2.0,
  "icFile":"$ic_file",
  "initialDt":1e-06,
  "nN":$nN,
  "nZ":$nZ,
  "saveFolder":"$save_folder/",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":0.05,
  "horizontalBoundaryConditions":"periodic",
  "isCudaEnabled":true,
  "threadsPerBlock_x":16,
  "threadsPerBlock_y":32
}
EOF

python3 tools/make_initial_conditions.py --output $ic_file --perturb_vorticity --n_modes $nN --n_gridpoints $nZ --modes 1

echo "==================== Building program"
#make clean
make gpu -j4

echo "==================== Starting program"
{ /usr/bin/time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log

echo "==================== Comparing results"

cat << EOF > $save_folder/benchmark.txt
0   6.34039E-01   -2.58587E+03  -2.38721E+02
1   -4.88785E-02  -1.88840E+02  -1.51367E+00
2   1.56819E-02   -6.27235E+01  -2.54914E+00
3   2.00381E-03   1.24415E+02   -2.78941E-03
4   -6.50265E-03  -6.35621E+01  1.65531E-01
5   2.02125E-03   -3.10272E+01  -1.12995E-02
6   1.61460E-03   4.82811E+01   -9.07794E-03
7   -1.36748E-03  -5.09181E+00  5.49947E-03
8   -1.61790E-04  -2.06828E+01  -2.39669E-03
9   5.84908E-04   9.11651E+00   -9.96605E-04
10  -9.91659E-05  6.34053E+00   1.49637E-03
11  -2.03341E-04  -5.56783E+00  -9.26890E-05
12  8.55685E-05   -1.29531E+00  -5.34289E-04
13  6.01665E-05   2.55982E+00   1.61954E-04
14  -4.49620E-05  3.95612E-03   1.51440E-04
15  -1.43399E-05  -1.02090E+00  -8.62517E-05
16  1.96831E-05   1.80896E-01   -3.48393E-05
17  1.86305E-06   3.69142E-01   3.52323E-05
18  -7.67912E-06  -1.32714E-01  5.39074E-06
19  6.51380E-07   -1.20063E-01  -1.25471E-05
20  2.70688E-06   7.07326E-02   3.60620E-07
EOF

comparison_results=$(python3 tools/print_variables.py $save_folder/dump0005.dat --max_print_mode 20 --constants $save_folder/constants.json | column -t | diff - $save_folder/benchmark.txt)

if [ "$comparison_results" ]; then
  echo "Simulation returned different results!"
  exit -1
else
  echo "Results match!"
  exit 0
fi
