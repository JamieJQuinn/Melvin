#pragma once

#include <constants.hpp>
#include <sim.hpp>

class DoubleDiffusiveSimulation : public Sim {
  public:
    Variable xi;
    Variable dXidt;

    DoubleDiffusiveSimulation(const Constants &c_in);
    void save();
    void load(const std::string &icFile);
    void reinit();
    void computeLinearDerivatives();
    void addAdvectionApproximation();
    void updateVars(real f=1.0);
    void advanceDerivatives();
    void runLinearStep();
    void runNonLinearStep(real f=1.0);
    void computeNonLinearDerivatives();
};
