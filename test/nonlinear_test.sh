#!/usr/bin/env bash

save_folder="test/benchmark"
nN=51
nZ=101

mkdir -p $save_folder
rm -f $save_folder/*

cat << EOF > $save_folder/constants.json
{
  "Pr":0.5,
  "Ra":1e6,
  "aspectRatio":3,
  "icFile":"$save_folder/ICn1nZ101nN51",
  "initialDt":3e-6,
  "nN":$nN,
  "nZ":$nZ,
  "saveFolder":"$save_folder/",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":0.05
}
EOF

constants_file=$save_folder/constants.json
python3 tools/make_initial_conditions.py --output $save_folder/ICn1nZ101nN51 --n_modes $nN --n_gridpoints $nZ --modes 1

echo "==================== Building program"
#make clean
make -j4 release

echo "==================== Starting program"
{ /usr/bin/time build/exe --constants $constants_file ; } 2>&1 | tee $save_folder/log

echo "==================== Comparing results"

cat << EOF > $save_folder/benchmark.txt
0   5.25484E-01   0.00000E+00   0.00000E+00
1   3.72859E-02   2.46098E+03   2.21268E+02
2   -2.86759E-02  -2.29944E+01  -2.59228E+00
3   2.99074E-02   1.03658E+03   6.41343E+01
4   4.08538E-03   1.04916E+02   5.77251E-01
5   2.82203E-02   9.09072E+02   2.77226E+01
6   7.86114E-03   4.84673E+01   1.67864E-01
7   3.18905E-02   8.05677E+02   1.34141E+01
8   5.16103E-03   -2.96553E+01  -3.91881E-01
9   3.60247E-02   6.78501E+02   6.91936E+00
10  2.13138E-03   -8.15562E+01  -5.58843E-01
11  3.64325E-02   5.39018E+02   3.77499E+00
12  6.21721E-04   -8.31905E+01  -4.61587E-01
13  3.31278E-02   4.19935E+02   2.16932E+00
14  1.57951E-04   -6.39944E+01  -2.96284E-01
15  2.78743E-02   3.31865E+02   1.30294E+00
16  1.20889E-04   -4.53150E+01  -1.68803E-01
17  2.23769E-02   2.67582E+02   8.12420E-01
18  1.78549E-04   -3.47940E+01  -9.82104E-02
19  1.76286E-02   2.17643E+02   5.23489E-01
20  2.59607E-04   -3.04392E+01  -6.53589E-02
EOF

comparison_results=$(python3 tools/print_variables.py $save_folder/dump0005.dat --max_print_mode 20 --constants $save_folder/constants.json | column -t | diff - $save_folder/benchmark.txt)

if [ "$comparison_results" ]; then
  echo "Simulation returned different results!"
  exit -1
else
  echo "Results match!"
  exit 0
fi
