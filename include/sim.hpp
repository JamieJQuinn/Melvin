#pragma once

#include <vector>

#include <thomas_algorithm.hpp>
#include <constants.hpp>
#include <precision.hpp>
#include <variable.hpp>

class Sim {
  public:
    bool modifydt; // has dt been modified?
    real t; // current time
    real dt; // current timestep
    int saveNumber; // Current save file

    // Kinetic Energy tracker
    std::vector<real> kineticEnergies;

    Constants c;

    // Variable arrays
    Variable psi; // Stream function (Psi)
    Variable omg; // Vorticity (Omega)
    Variable tmp; // Temperature

    Variable dTmpdt; // d/dt of temperature
    Variable dOmgdt; // d/dt of vorticty

    Variable xi; // Double diffusive constituent
    Variable dXidt; // d/dt of xi

    std::vector<Variable*> variableList; // for operations like saving, loading and initialising

    ThomasAlgorithm *thomasAlgorithm;

    Sim(const Constants &c_in);
    ~Sim();

    // Helper functions
    void initialiseThomasAlgorithm();
    void printMaxOf(real *a, std::string name) const;
    void printBenchmarkData() const;
    void save();
    std::string createSaveFilename();
    void load(const std::string &icFile);
    void reinit();
    void calcKineticEnergy();
    real calcKineticEnergyForMode(int n);
    void saveKineticEnergy();
    real isFinished();

    // Derivative calculations
    void computeLinearDerivatives();
    void computeLinearTemperatureDerivative();
    void computeLinearVorticityDerivative();
    void computeLinearXiDerivative();
    void addAdvectionApproximation();

    void computeNonlinearDerivatives();
    void computeNonlinearDerivative(Variable &dVardt, const Variable &var);
    void computeNonlinearTemperatureDerivative();
    void computeNonlinearXiDerivative();
    void computeNonlinearVorticityDerivative();

    // Simulation functions
    void solveForPsi();
    void applyBoundaryConditions();
    void updateVars(real f=1.0);
    void advanceDerivatives();

    void runNonLinear();
    void runNonLinearStep(real f=1.0);

    // Linear critical Rayleigh functions
    int testCriticalRayleigh();
    bool isCritical(int nCrit);
    real findCriticalRa(int nCrit);
    void runLinearStep();
};
