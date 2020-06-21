#pragma once

#include <vector>

#include <thomas_algorithm.hpp>
#include <constants.hpp>
#include <precision.hpp>
#include <variable.hpp>
#include <variables.hpp>
#include <kinetic_energy_tracker.hpp>

#ifdef CUDA
#include <variable_gpu.hpp>
#include <thomas_algorithm_gpu.hpp>
#endif

class Sim {
  public:
    bool modifydt; // has dt been modified?
    real t; // current time
    real dt; // current timestep

    // Kinetic Energy tracker
    KineticEnergyTracker keTracker;

    Constants c;

    // Variable arrays
    Variables<Variable> vars;
    Variable nonlinearSineTerm, nonlinearCosineTerm;

    ThomasAlgorithm *thomasAlgorithm;

    Sim(const Constants &c_in);
    ~Sim();

    // Helper functions
    void printMaxOf(real *a, std::string name) const;
    void printBenchmarkData() const;
    real isFinished();

    // Derivative calculations
    void computeLinearDerivatives();
    void computeLinearTemperatureDerivative();
    void computeLinearVorticityDerivative();
    void computeLinearXiDerivative();
    void addAdvectionApproximation();

    void computeNonlinearDerivatives();
    void applyPhysicalBoundaryConditions();
    void computeNonlinearDerivative(Variable &dVardt, const Variable &var);
    void computeNonlinearTemperatureDerivative();
    void computeNonlinearXiDerivative();
    void computeNonlinearVorticityDerivative();

    // Simulation functions
    void solveForPsi();
    void applyBoundaryConditions();
    void applyTemperatureBoundaryConditions();
    void applyVorticityBoundaryConditions();
    void applyPsiBoundaryConditions();
    void applyXiBoundaryConditions();

    void runNonLinear();
    void runNonLinearStep(real f=1.0);

    // Linear critical Rayleigh functions
    int testCriticalRayleigh();
    bool isCritical(int nCrit);
    real findCriticalRa(int nCrit);
    void runLinearStep();
};
