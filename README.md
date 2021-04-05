# Melvin

Melvin is a 2D, pseudo-spectral, adaptive-time implementation of
* Rayleigh-Benard convection (RBC)
* Double-diffusive convection (DDC)

The code has been parallelised on the CPU via OpenMP and FFTW and on the GPU via CUDA and cuFFT. Not all features and boundary conditions are currently implemented.

## Why is it called Melvin?

[Melvin Stern](https://en.wikipedia.org/wiki/Melvin_Stern) was a pioneer in oceanographic fluid dynamics, particularly in the study of salt-fingering. This code is a personal monument to a great scientist.

## Dependencies

### Running on CPU

- FFTW 3
- gcc with OpenMP support (tested with gcc 6 and 9)

### Running on GPU

Identical dependencies as CPU version plus:

- CUDA development and runtime libraries

## Usage

The various run scripts available in `run_scripts` will automatically build and run the simulation. They are currently set up to run a number of different cases including

- Rayleigh-Benard convection
- Double-diffusive salt-fingering with and without periodic boundary conditions
- Rayleigh-Taylor instability 

These scripts may be used to extend the simulations to other parameter choices, boundary conditions and initial conditions.

For manual building, running `make` will automatically build the OpenMP enabled CPU version of the code. `make profile` will include profiling information in the binary while `make debug` will enable assertions and include symbols in the binary.

To compile for GPU run `make gpu`. There is a corresponding `make gpu-debug` for debuggin purposes.

## Tests

### Regression tests

Currently the regression tests available are:

- Linear:
  - Standard Rayleigh-Benard convection
  - Salt-fingering
- Nonlinear:
  - Standard Rayleigh-Benard convection
  - Combined compositional + thermal convection
  - Salt-fingering

The linear tests will run two linear simulations either side of the analytically calculated critical Rayleigh number.

The nonlinear RBC test will run to completion and compare the output to the known test found in Glatzmaier, page 46.

All tests can be found in the `test` folder.

### Unit tests

The unit-testing framework currently used is the header-only library [Catch](https://github.com/catchorg/Catch2). 

To run the unit tests, run `make test`. The GPU tests can be run with `make gpu-test`.

## Boundary conditions

For **RBC only** there is periodic and impermeable boundary conditions available in either direction. Due to constraints on the FFT algorithm, the periodic boundary conditions are faster for most resolutions.

## Performance

Due to the nature of the fast-Fourier transform algorithm, the running performance is best when the horizontal real-space resolution is a factor of small primes. Powers of two are preferrable.

## References

> Glatzmaier: Introduction to Modeling Convection in Stars and Planets; Gary A. Glatzmaier; 2014
