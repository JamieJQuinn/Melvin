#pragma once

#include <variable.hpp>

class VariableGPU: public Variable {
  public:
    void initialiseData(real initialValue = 0.0);
    VariableGPU(const Constants &c_in, int totalSteps_in=1);
    ~VariableGPU();
};
