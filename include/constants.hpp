#pragma once

#include <string>

class Constants {
  public:
    int nZ;
    int nN;
    double initialDt;
    double Ra;
#ifdef DDC
    double RaXi;
    double tau;
#endif
    double Pr;
    int aspectRatio;
    double timeBetweenSaves;
    double totalTime;

    // Derived ants
    double dz;
    double dx;
    int nX;
    double oodz2;

    // I/O files
    std::string saveFolder;
    std::string icFile; // initial conditions

    Constants();
    Constants(const std::string &input);
    Constants(int nZ_in, int nN_in, double initialDt_in, double Ra_in, double Pr_in, int aspectRatio_in, double timeBetweenSaves_in, double totalTime_in, std::string &saveFolder_in, std::string &icFile_in);

    void calculateDerivedConstants();

    void print() const;
    bool isValid() const;

    void readJson(const std::string &filePath);
    void writeJson(const std::string &filePath) const;
};
