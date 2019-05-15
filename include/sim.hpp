#pragma once

#include <thomas_algorithm.hpp>
#include <constants.hpp>
#include <precision.hpp>

class Sim {
  public:
    bool modifydt; // has dt been modified?
    int current; // which array is the current derivative in?
    real t; // current time
    real dt; // current timestep
    int saveNumber; // Current save file
    int KEsaveNumber; // Current kinetic energy save file

    // Gradients HACK - this needs refactored somewhere - into Constants?
    real tmpGrad;
    real xiGrad;

    // Kinetic Energy tracker
    real kePrev;
    real keCurrent;

    const Constants c;

    // Variable arrays
    real * psi; // Stream function (Psi)
    real * omg; // Vorticity (Omega)
    real * tmp; // Temperature
#ifdef DDC
    real * xi;  // Salt concentration
    real * dXidt; // d/dt of salt concentration
#endif

    real * dOmgdt; // d/dt of vorticty
    real * dTmpdt; // d/dt of temperature

    ThomasAlgorithm *thomasAlgorithm;

    Sim(const Constants &c_in);
    ~Sim();

    // Helper functions
    void printMaxOf(real *a, std::string name);
    void printBenchmarkData();
    void save();
    void load(real* tmp, real* omg, real* psi, const std::string &icFile);
    void reinit();
    real calcKineticEnergy();
    real calcKineticEnergyForMode(int n);
    void saveKineticEnergy();

    // Simulation functions
    void updateTmpAndOmg(real f);
#ifdef DDC
    void updateXi(real f);
#endif
    void computeLinearDerivatives(int linearSim = 1);
    void computeNonLinearDerivatives();
    void solveForPsi();

    // Runs the linear simulation
    real runLinear(int);
    void runNonLinear();
};
