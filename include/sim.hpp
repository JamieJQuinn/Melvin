#pragma once

#include "thomas_algorithm.hpp"
#include "constants.hpp"

class Sim {
  public:
    bool modifydt; // has dt been modified?
    int current; // which array is the current derivative in?
    double t; // current time
    double dt; // current timestep
    int saveNumber; // Current save file
    int KEsaveNumber; // Current kinetic energy save file

    // Gradients HACK - this needs refactored somewhere - into Constants?
    double tmpGrad;
    double xiGrad;

    // Kinetic Energy tracker
    double kePrev;
    double keCurrent;

    const Constants c;

    // Variable arrays
    double * psi; // Stream function (Psi)
    double * omg; // Vorticity (Omega)
    double * tmp; // Temperature
#ifdef DDC
    double * xi;  // Salt concentration
    double * dXidt; // d/dt of salt concentration
#endif

    double * dOmgdt; // d/dt of vorticty
    double * dTmpdt; // d/dt of temperature

    ThomasAlgorithm *thomasAlgorithm;

    Sim(const Constants &c_in);
    ~Sim();

    // Helper functions
    void printMaxOf(double *a, std::string name);
    void printBenchmarkData();
    void save();
    void load(double* tmp, double* omg, double* psi, const std::string &icFile);
    void reinit();
    double calcKineticEnergy();
    double calcKineticEnergyForMode(int n);
    void saveKineticEnergy();

    // Simulation functions
    void updateTmpAndOmg(double f);
#ifdef DDC
    void updateXi(double f);
#endif
    void computeLinearDerivatives(int linearSim = 1);
    void computeNonLinearDerivatives();
    void solveForPsi();

    // Runs the linear simulation
    double runLinear(int);
    void runNonLinear();
};
