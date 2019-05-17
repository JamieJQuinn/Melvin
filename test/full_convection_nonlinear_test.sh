#!/usr/bin/env bash

save_folder="test/benchmark"

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.js
{
  "Pr":0.5,
  "Ra":50000,
  "RaXi":100000,
  "tau":0.01,
  "aspectRatio":3,
  "icFile":"$save_folder/ICn1nZ101nN51_SF",
  "initialDt":3e-6,
  "nN":51,
  "nZ":101,
  "saveFolder":"$save_folder/",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":true,
  "totalTime":0.05
}
EOF

constants_file=$save_folder/constants.js
python tools/make_initial_conditions.py --output $save_folder/ICn1nZ101nN51_SF --n_modes 51 --n_gridpoints 101 --modes 1 --combined_convection

echo "==================== Building program"
make clean
make release

echo "==================== Starting program"
{ time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log

echo "==================== Comparing results"

cat << EOF > $save_folder/benchmark.txt
0   6.80000E-01   0.00000E+00   0.00000E+00
1   -6.52901E-04  -4.35834E-01  -3.97463E-02
2   2.18389E-13   2.41244E-10   4.56874E-12
3   -1.52664E-16  -1.50111E-13  -9.54637E-15
4   -3.52353E-20  -1.21786E-16  -1.75760E-18
5   2.79230E-23   5.65021E-20   1.89796E-21
6   1.57642E-27   8.70101E-24   7.02881E-26
7   -1.27449E-30  -1.54184E-27  -5.88412E-29
8   -5.14961E-35  -2.80953E-30  -2.72398E-32
9   -1.06029E-39  1.63078E-33   1.85969E-35
10  9.66997E-41   -9.83341E-37  -7.53473E-39
11  -8.30608E-44  7.29499E-40   5.18361E-42
12  7.32243E-47   -1.26670E-43  -3.18862E-45
13  -5.99771E-50  -1.53720E-46  2.02327E-48
14  4.95632E-53   -1.69041E-48  -1.42706E-51
15  -3.90770E-56  2.81775E-51   1.67208E-54
16  2.94733E-59   3.92425E-54   -1.14833E-57
17  -2.33392E-62  -8.07590E-57  -1.40195E-60
18  2.43816E-65   -1.27397E-59  2.21301E-63
19  -2.47793E-68  1.94826E-62   3.21177E-66
#20  7.38886E-72   5.11045E-65   -6.49594E-69
EOF

comparison_results=$(python3 tools/print_variables.py $save_folder/dump0005.dat --max_print_mode 20 --n_modes 51 --n_gridpoints 101 | column -t | diff - $save_folder/benchmark.txt)

if [ "$comparison_results" ]; then
  echo "Simulation returned different results!"
  exit -1
else
  echo "Results match!"
  exit 0
fi
