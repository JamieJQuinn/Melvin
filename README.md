# Introduction-To-Modelling

This project includes implementations of
* Linear Rayleigh-Benard convection
* Nonlinear Rayleigh-Benard convection
* Linear double-diffusive convection

## Integration Test

Running 
```
test/integration_test.sh
```
will
1. build the nonlinear simulation 
2. run (inputs may be found in script)
3. compare the output to the "correct" output found in `data/benchmark_t5.txt`

This output is correct in so far as running to T=3 will compare exactly with the solutions found in Glatzmaier[^1], page 46.

[^1]: Introduction to Modeling Convection in Stars and Planets; Gary A. Glatzmaier; 2014
