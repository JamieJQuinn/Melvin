#pragma once

#include <string>
#include <precision.hpp>

class Constants {
  public:
    int nZ;
    int nN;
    real initialDt;
    real Ra;
#ifdef DDC
    real RaXi;
    real tau;
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
    Constants(int nZ_in, int nN_in, real initialDt_in, real Ra_in, real Pr_in, int aspectRatio_in, real timeBetweenSaves_in, real totalTime_in, std::string &saveFolder_in, std::string &icFile_in);

    void calculateDerivedConstants();

    void print() const;
    bool isValid() const;

    void readJson(const std::string &filePath);
    void writeJson(const std::string &filePath) const;
};
