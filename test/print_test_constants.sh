#!/usr/bin/env bash

cat << EOF > test_constants.json
{
  "Pr":0.5,
  "Ra":1000000,
  "aspectRatio":3,
  "icFile":"ICn1nZ101nN51",
  "initialDt":3e-06,
  "nN":64,
  "nZ":128,
  "saveFolder":".",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":0.05
}
EOF

cat << EOF > test_constants_gpu.json
{
  "Pr":0.5,
  "Ra":1e4,
  "RaXi":1e5,
  "tau":0.01,
  "aspectRatio":3,
  "icFile":"ICn1nZ101nN51",
  "initialDt":3e-06,
  "nN":64,
  "nZ":128,
  "saveFolder":".",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":true,
  "isCudaEnabled":true,
  "threadsPerBlock_x":16,
  "threadsPerBlock_y":32,
  "totalTime":0.05
}
EOF
