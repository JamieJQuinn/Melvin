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

    void computeNonlinearDerivativeSpectralTransform(VariableGPU &dVardt, const VariableGPU &var);
    void computeNonlinearDerivatives();
    void computeNonlinearTemperatureDerivative();
    void computeNonlinearXiDerivative();
    void computeNonlinearVorticityDerivative();

    void solveForPsi();

    void addAdvectionApproximation();

    void runLinearStep();
    void runNonLinearStep(real f=1.0);

    void runNonLinear();

    VariableGPU nonlinearTerm;
  private:
    ThomasAlgorithmGPU *thomasAlgorithm;
};
