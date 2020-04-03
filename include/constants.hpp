#pragma once

#include <string>
#include <precision.hpp>
#include <boundary_conditions.hpp>

class Constants {
  public:
    int nZ;
    int nN;
    int nX;
    int nG;
    real initialDt;
    real Ra;
    real RaXi;
    real tau;
    real Pr;
    real aspectRatio;
    real timeBetweenSaves;
    real totalTime;
    bool isNonlinear;
    bool isDoubleDiffusion;
    bool isCudaEnabled;

    bool isPhysicalResSpecfified;

    // Boundary conditions
    BoundaryConditions verticalBoundaryConditions;
    BoundaryConditions horizontalBoundaryConditions;
    real temperatureGradient;
    real salinityGradient;

    // cuda stuff
    int threadsPerBlock_x;
    int threadsPerBlock_y;

    // Derived ants
    real dz;
    real dx;
    real oodz2;
    real oodz;
    real oodx;
    mode wavelength;

    // I/O files
    std::string saveFolder;
    std::string icFile; // initial conditions

    Constants();
    ~Constants();
    Constants(const std::string &input);

    void calculateDerivedConstants();

    void print() const;
    bool isValid() const;

    void readJson(const std::string &filePath);
    void writeJson(const std::string &filePath) const;

  private:
    std::string verticalBoundaryConditions_in;
    std::string horizontalBoundaryConditions_in;
};
