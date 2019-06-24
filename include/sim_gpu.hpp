#pragma once

#include <variables.hpp>
#include <variable_gpu.hpp>
#include <thomas_algorithm_gpu.hpp>
#include <kinetic_energy_tracker.hpp>
#include <constants.hpp>

class SimGPU {
  public:
    bool modifydt; // has dt been modified?
    real t; // current time
    real dt; // current timestep

    KineticEnergyTracker keTracker;

    Constants c;

    // Variable arrays
    Variables<VariableGPU> vars;

    SimGPU(const Constants &c_in);
    ~SimGPU();

    void computeLinearDerivatives();
    void computeLinearTemperatureDerivative();
    void computeLinearVorticityDerivative();
    void computeLinearXiDerivative();

    void solveForPsi();

    void addAdvectionApproximation();

    void runLinearStep();

  private:
    ThomasAlgorithmGPU *thomasAlgorithm;
};
