#pragma once

#include <variable.hpp>
#include <precision.hpp>

class VariableGPU: public Variable {
  const int threadsPerBlock_x;
  const int threadsPerBlock_y;

  public:
    void initialiseData(real initialValue = 0.0);
    VariableGPU(const Constants &c_in, int totalSteps_in=1, bool useSinTransform_in = true);
    void update(const VariableGPU& dVardt, const real dt, const real f=1.0);

    gpu_mode * getCurrent();
    const gpu_mode * getCurrent() const;
    gpu_mode* getPrevious();
    const gpu_mode* getPrevious() const;
    gpu_mode * getPlus(int nSteps=0);
    const gpu_mode* getPlus(int nSteps) const;

    void readFromFile(std::ifstream& file);
    ~VariableGPU();
};
