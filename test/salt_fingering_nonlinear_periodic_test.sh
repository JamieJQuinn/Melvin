#!/usr/bin/env bash

save_folder="test/benchmark"

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.json
{
  "Pr":1,
  "Ra":1.1e4,
  "RaXi":1e6,
  "tau":1e-2,
  "aspectRatio":1.41421356237,
  "initialDt":3e-6,

  "nN":85,
  "nZ":101,
  "nX":256,

  "horizontalBoundaryConditions": "periodic",
  "verticalBoundaryConditions": "periodic",
  "temperatureGradient":1,
  "salinityGradient":1,

  "icFile":"$save_folder/ICn1nZ101nN51_SF",
  "saveFolder":"$save_folder/",

  "timeBetweenSaves":0.1,
  "totalTime":0.5,

  "isNonlinear":true,
  "isDoubleDiffusion":false
}
EOF



constants_file=$save_folder/constants.json
python tools/make_initial_conditions.py --output $save_folder/ICn1nZ101nN51_SF --n_modes 85 --n_gridpoints 101 --modes 2 --amp 0.001 --salt_fingering --periodic

#--modes 1 25 --salt_fingering --periodic

echo "==================== Building program"
#make clean
make release

echo "==================== Starting program"
{ /usr/bin/time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log

echo "==================== Comparing results"
cat << EOF > $save_folder/benchmark.txt
0   3.20258E-01   0.00000E+00   0.00000E+00   3.20811E-01
1   1.67816E-02   -7.18441E+00  -3.63967E-01  2.96803E-02
2   -1.11473E-07  1.97964E-04   2.18067E-06   -9.05513E-07
3   1.10592E-09   -3.02279E-06  -3.62231E-08  5.82295E-09
4   7.34376E-12   -6.25388E-08  -2.84496E-10  1.77342E-10
5   -9.01203E-14  1.02746E-09   4.20199E-12   -1.93706E-12
6   -7.84467E-16  7.63426E-12   1.74609E-14   -2.93554E-14
7   -2.46271E-17  3.31325E-13   6.71455E-16   -7.90845E-16
8   -1.31952E-18  3.01888E-14   4.22181E-17   -8.53801E-17
9   -6.13177E-20  1.14761E-15   1.54404E-18   -2.38332E-18
10  -2.25607E-21  6.68484E-17   6.58812E-20   -1.43916E-19
11  3.46374E-21   -2.62098E-17  -4.62705E-20  3.75951E-19
12  8.42183E-20   4.04984E-16   -8.77365E-19  -9.28981E-18
13  1.81032E-18   1.05262E-14   -1.82197E-17  -3.41013E-16
14  3.61217E-17   -4.80345E-13  -4.83861E-16  -3.65510E-16
15  6.28952E-16   -2.13152E-11  -1.01927E-14  1.46022E-13
16  9.16184E-15   -2.66591E-10  -1.28917E-13  2.35037E-12
17  1.20042E-13   1.66681E-09   -8.11035E-13  -1.51263E-11
18  1.72746E-12   7.18115E-08   -8.88962E-12  -9.52040E-10
19  2.82705E-11   -1.93457E-07  -4.50066E-10  -1.07569E-08
20  4.19630E-10   -3.51607E-05  -1.19513E-08  4.76471E-08
EOF

comparison_results=$(python3 tools/print_variables.py $save_folder/dump0005.dat --max_print_mode 20 --constants $save_folder/constants.json | column -t | diff - $save_folder/benchmark.txt)

if [ "$comparison_results" ]; then
  echo "Simulation returned different results!"
  exit -1
else
  echo "Results match!"
  exit 0
fi
