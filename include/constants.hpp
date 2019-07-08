#pragma once

#include <string>
#include <precision.hpp>

class Constants {
  public:
    int nZ;
    int nN;
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

    // cuda stuff
    int threadsPerBlock_x;
    int threadsPerBlock_y;

    // Derived ants
    real dz;
    real dx;
    int nX;
    real oodz2;
    real oodz;

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
};
