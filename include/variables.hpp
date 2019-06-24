#pragma once

#include <vector>

#include <variable.hpp>

class Variables {
  public:
    Variable tmp;
    Variable xi;
    Variable omg;
    Variable psi;
    Variable dTmpdt;
    Variable dOmgdt;
    Variable dXidt;

    int saveNumber;

    const Constants c;

    Variables(const Constants &c_in);

    std::vector<Variable*> variableList; // for operations like saving, loading and initialising

    void save();
    void load(const std::string &icFile);
    void reinit(const real value = 0.0);

    void updateVars(const real dt, const real f = 1.0);
    void advanceDerivatives();
  private:
    std::string createSaveFilename();
};
