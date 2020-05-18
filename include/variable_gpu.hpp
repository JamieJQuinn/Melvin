#pragma once

#include <variable.hpp>
#include <precision.hpp>

__device__ int calcIndex(int n, int k);

class VariableGPU: public Variable {
  const int threadsPerBlock_x;
  const int threadsPerBlock_y;

  gpu_mode * data_d;

  public:
    void initialiseData(mode initialValue = 0.0);
    VariableGPU(const Constants &c_in, int totalSteps_in=1, bool useSinTransform_in = true);
    void update(const VariableGPU& dVardt, const real dt, const real f=1.0);
    void fill(const mode value);

    void copyToDevice();
    void copyToHost();

    gpu_mode * getCurrent();
    const gpu_mode * getCurrent() const;
    gpu_mode* getPrevious();
    const gpu_mode* getPrevious() const;
    gpu_mode * getPlus(int nSteps=0);
    const gpu_mode* getPlus(int nSteps) const;

    void readFromFile(std::ifstream& file);
    void writeToFile(std::ofstream& file);
    ~VariableGPU();
};
