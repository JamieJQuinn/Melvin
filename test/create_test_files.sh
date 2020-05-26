#!/usr/bin/env bash

n_modes=170
n_transform_points=512
n_gridpoints=256

cat << EOF > test_constants.json
{
  "Pr":0.5,
  "Ra":1000000,
  "aspectRatio":3,
  "icFile":"initial_conditions.dat",
  "initialDt":3e-06,
  "nN":${n_modes},
  "nZ":${n_gridpoints},
  "saveFolder":"./",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":0.05,

  "isCudaEnabled":true,
  "threadsPerBlock_x":16,
  "threadsPerBlock_y":32
}
EOF

cat << EOF > test_constants_cpu.json
{
  "Pr":0.5,
  "Ra":1000000,
  "aspectRatio":3,
  "icFile":"initial_conditions.dat",
  "initialDt":3e-06,
  "nN":${n_modes},
  "nZ":${n_gridpoints},
  "saveFolder":"./",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":0.05,

  "isCudaEnabled":false,
  "threadsPerBlock_x":16,
  "threadsPerBlock_y":32
}
EOF

cat << EOF > test_constants_periodic_gpu.json
{
  "Pr":0.5,
  "Ra":1000000,
  "aspectRatio":3,
  "icFile":"initial_conditions.dat",
  "initialDt":3e-06,
  "nN":${n_modes},
  "nX":${n_transform_points},
  "nZ":${n_gridpoints},
  "saveFolder":"./",
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

cat << EOF > test_constants_periodic.json
{
  "Pr":0.5,
  "Ra":1000000,
  "aspectRatio":3,
  "icFile":"initial_conditions.dat",
  "initialDt":3e-06,
  "nN":${n_modes},
  "nX":${n_transform_points},
  "nZ":${n_gridpoints},
  "saveFolder":"./",
  "timeBetweenSaves":0.01,
  "isNonlinear":true,
  "isDoubleDiffusion":false,
  "totalTime":0.05,
  "horizontalBoundaryConditions":"periodic"
}
EOF

cat << EOF > test_constants_ddc_gpu.json
{
  "Pr":1,
  "Ra":1e6,
  "RaXi":1e7,
  "tau":1e-2,
  "aspectRatio":1.41421356237,
  "initialDt":3e-6,

  "nN":${n_modes},
  "nZ":${n_gridpoints},

  "icFile":"initial_conditions_ddc.dat",
  "saveFolder":"./",

  "timeBetweenSaves":0.01,
  "totalTime":0.05,

  "isNonlinear":true,
  "isDoubleDiffusion":true,

  "isCudaEnabled":true,
  "threadsPerBlock_x":16,
  "threadsPerBlock_y":32
}
EOF

	python3 ../tools/make_initial_conditions.py --output initial_conditions_ddc.dat --salt_fingering --n_modes ${n_modes} --n_gridpoints ${n_gridpoints} --modes $(seq 1 63)
	python3 ../tools/make_initial_conditions.py --output initial_conditions.dat --n_modes ${n_modes} --n_gridpoints ${n_gridpoints} --modes $(seq 1 63)

