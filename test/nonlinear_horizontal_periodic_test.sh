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
  "totalTime":0.5,
  "horizontalBoundaryConditions":"periodic"
}
EOF

python3 tools/make_initial_conditions.py --output $save_folder/initial_conditions.dat --n_modes $nN --n_gridpoints $nZ --modes 1

echo "==================== Building program"
#make clean
make -j4 release

echo "==================== Starting program"
{ /usr/bin/time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log

echo "==================== Comparing results"

cat << EOF > $save_folder/benchmark.txt
0   5.01752E-01   5.15469E-03   5.22322E-04
1   3.43080E-02   3.48057E+03   3.10385E+02
2   -2.17143E-03  1.60229E+01   -1.06352E+00
3   2.96567E-02   2.01443E+03   6.89928E+01
4   1.67522E-03   -5.28176E-01  -6.42442E-01
5   3.18523E-02   1.50814E+03   2.21529E+01
6   1.13490E-03   -3.69291E+01  -5.61085E-01
7   3.28215E-02   1.06826E+03   8.33806E+00
8   1.19618E-03   -5.24982E+01  -4.02572E-01
9   3.18270E-02   7.83373E+02   3.77375E+00
10  1.64479E-03   -5.97916E+01  -2.66261E-01
11  2.95327E-02   6.07688E+02   1.98005E+00
12  2.06124E-03   -6.35041E+01  -1.83028E-01
13  2.65150E-02   4.85509E+02   1.13809E+00
14  2.31359E-03   -6.27857E+01  -1.28777E-01
15  2.31302E-02   3.91675E+02   6.91913E-01
16  2.40111E-03   -5.72342E+01  -8.86911E-02
17  1.96574E-02   3.16954E+02   4.37192E-01
18  2.35614E-03   -4.85426E+01  -5.90670E-02
19  1.63420E-02   2.56931E+02   2.84354E-01
20  2.22164E-03   -3.88901E+01  -3.81715E-02
EOF

comparison_results=$(python3 tools/print_variables.py $save_folder/dump0005.dat --max_print_mode 20 --constants $save_folder/constants.json | column -t | diff - $save_folder/benchmark.txt)

if [ "$comparison_results" ]; then
  echo "Simulation returned different results!"
  exit -1
else
  echo "Results match!"
  exit 0
fi
