#!/usr/bin/env bash

save_folder="test/benchmark"

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
  "totalTime":0.05
}
EOF

constants_file=$save_folder/constants.js
python3 tools/make_initial_conditions.py --output $save_folder/ICn1nZ101nN51 --n_modes 51 --n_gridpoints 101 --modes 1

echo "==================== Building program"
#make clean
make release

echo "==================== Starting program"
{ /usr/bin/time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log

echo "==================== Comparing results"

cat << EOF > $save_folder/benchmark.txt
0   5.25632E-01   0.00000E+00   0.00000E+00
1   3.72935E-02   2.45535E+03   2.20950E+02
2   -2.90106E-02  -2.22633E+01  -2.57855E+00
3   2.98495E-02   1.03644E+03   6.40933E+01
4   3.83443E-03   1.05693E+02   6.29383E-01
5   2.83494E-02   9.07609E+02   2.77101E+01
6   7.76393E-03   4.63140E+01   1.95161E-01
7   3.19169E-02   8.03862E+02   1.34052E+01
8   5.13458E-03   -3.15273E+01  -3.92410E-01
9   3.60185E-02   6.75509E+02   6.90593E+00
10  2.09485E-03   -8.36919E+01  -5.68721E-01
11  3.64468E-02   5.35579E+02   3.76424E+00
12  5.52325E-04   -8.56798E+01  -4.74691E-01
13  3.31556E-02   4.17356E+02   2.16357E+00
14  6.01454E-05   -6.76310E+01  -3.10357E-01
15  2.79294E-02   3.30956E+02   1.30178E+00
16  6.03942E-06   -5.02731E+01  -1.82265E-01
17  2.24698E-02   2.68417E+02   8.14412E-01
18  6.61424E-05   -4.04443E+01  -1.09780E-01
19  1.77579E-02   2.19639E+02   5.27082E-01
20  1.65159E-04   -3.56803E+01  -7.44043E-02
EOF

comparison_results=$(python3 tools/print_variables.py $save_folder/dump0005.dat --max_print_mode 20 --n_modes 51 --n_gridpoints 101 | column -t | diff - $save_folder/benchmark.txt)

if [ "$comparison_results" ]; then
  echo "Simulation returned different results!"
  exit -1
else
  echo "Results match!"
  exit 0
fi
