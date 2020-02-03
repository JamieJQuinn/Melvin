#!/usr/bin/env bash

cat << EOF > test_constants.json
{
  "Pr":0.5,
  "Ra":1000000,
  "aspectRatio":3,
  "icFile":"ICn1nZ128nN64",
  "initialDt":3e-06,
  "nN":64,
  "nZ":128,
  "saveFolder":".",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":0.05,

  "isCudaEnabled":true,
  "threadsPerBlock_x":16,
  "threadsPerBlock_y":32
}
EOF

cat << EOF > test_constants_ddc.json
{
  "Pr":1,
  "Ra":1e6,
  "RaXi":1e7,
  "tau":1e-2,
  "aspectRatio":1.41421356237,
  "initialDt":3e-6,

  "nN":64,
  "nZ":128,

  "icFile":"ICn1nZ128nN64_SF",
  "saveFolder":".",

  "timeBetweenSaves":0.01,
  "totalTime":0.05,

  "isNonlinear":true,
  "isDoubleDiffusion":true,

  "isCudaEnabled":true,
  "threadsPerBlock_x":16,
  "threadsPerBlock_y":32
}
EOF
