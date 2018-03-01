#pragma once

class Sim {
  public:
    // Defined constants
    int nZ;
    int nN;
    double dt;
    double Ra;
#ifdef DDC
    double RaXi;
    double tau;
#endif
    double Pr;
          int a;
    double timeBetweenSaves;
    bool modifydt;
          int current;
    double t;
    double totalTime;
    double tmpGrad;
#ifdef DDC
    double xiGrad;
#endif

    // Derived constants
    double dz;
    double dx;
    int nX;
    double oodz2;
    int saveNumber;
    int KEsaveNumber;
    double timeBetweenKESaves;

    // Kinetic Energy tracker
    double kePrev;
    double keCurrent;

    // Save Folder
    std::string saveFolder;

    // Initial condition file
    std::string icFile;

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

    double * wk1; // Used in Thomas Algorithm
    double * wk2;
    double * sub;

    // Constructor
    Sim(int nZ, int nN, double dt,
    double Ra, double Pr, int a,
#ifdef DDC
    double RaXi, double tau,
#endif
    double timeBetweenSaves, bool modifydt,
          int current, double t, double totalTime,
    std::string saveFolder, std::string icFile);

    void init(int nZ, int nN, double dt,
    double Ra, double Pr, int a,
#ifdef DDC
    double RaXi, double tau,
#endif
    double timeBetweenSaves, bool modifydt,
          int current, double t, double totalTime,
    std::string saveFolder, std::string icFile);
    // Destructor
    ~Sim();

    // Helper functions
    void triDiagonalSolver(const int nZ,
             const double *rhs, double *sol, const double *sub,
             const double * wk1, const double *wk2);
    void formTridiArrays ( const int nZ,
      const double *sub, const double *dia, const double *sup,
      double * wk1, double *wk2);
    void printMaxOf(double *a, std::string name);
    void printBenchmarkData();
    void save();
    void load(double* tmp, double* omg, double* psi, std::string icFile);
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
