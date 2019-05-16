#pragma once

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
    int KEsaveNumber; // Current kinetic energy save file

    // Kinetic Energy tracker
    real kePrev;
    real keCurrent;

    Constants c;

    // Variable arrays
    Variable psi; // Stream function (Psi)
    Variable omg; // Vorticity (Omega)
    Variable tmp; // Temperature

    Variable dTmpdt; // d/dt of temperature
    Variable dOmgdt; // d/dt of vorticty

    ThomasAlgorithm *thomasAlgorithm;

    Sim(const Constants &c_in);
    ~Sim();

    // Helper functions
    void printMaxOf(real *a, std::string name) const;
    void printBenchmarkData() const;
    void save();
    void load(const std::string &icFile);
    void reinit();
    real calcKineticEnergy();
    real calcKineticEnergyForMode(int n);
    void saveKineticEnergy();
    real isFinished();

    // Simulation functions
    void computeLinearDerivatives();
    void addAdvectionApproximation();
    void computeNonLinearDerivatives();
    void solveForPsi();
    void applyBoundaryConditions();
    void updateVars(real f=1.0);
    void advanceDerivatives();

    void runNonLinear();

    // Linear critical Rayleigh functions
    void initialLinearConditions();
    bool isCritical(int nCrit);
    real findCriticalRa(int nCrit);
    virtual void runLinearStep();
};
