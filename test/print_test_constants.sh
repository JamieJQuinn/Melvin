#!/usr/bin/env bash

cat << EOF > test_constants.json
{
  "Pr":0.5,
  "Ra":1000000,
  "aspectRatio":3,
  "icFile":"ICn1nZ101nN51",
  "initialDt":3e-06,
  "nN":5,
  "nZ":10,
  "saveFolder":".",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":0.05
}
EOF
