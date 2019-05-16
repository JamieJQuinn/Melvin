# Introduction-To-Modelling

This project includes implementations of
* Linear and nonlinear Rayleigh-Benard convection
* Linear and nonlinear salt-fingering

## Tests

Currently the tests available are:

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

## References

> Glatzmaier: Introduction to Modeling Convection in Stars and Planets; Gary A. Glatzmaier; 2014
