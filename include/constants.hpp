#pragma once

#include <string>
#include <precision.hpp>

class Constants {
  public:
    int nZ;
    int nN;
    real initialDt;
    real Ra;
    int tempGrad;
#ifdef DDC
    real RaXi;
    real tau;
    int xiGrad;
#endif
    real Pr;
    int aspectRatio;
    real timeBetweenSaves;
    real totalTime;

    // Derived ants
    real dz;
    real dx;
    int nX;
    real oodz2;

    // I/O files
    std::string saveFolder;
    std::string icFile; // initial conditions

    Constants();
    Constants(const std::string &input);

    void calculateDerivedConstants();

    void print() const;
    bool isValid() const;

    void readJson(const std::string &filePath);
    void writeJson(const std::string &filePath) const;
};
