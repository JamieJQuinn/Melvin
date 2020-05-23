#pragma once

#include <variable.hpp>
#include <precision.hpp>

#include <cufft.h>

__device__ int calcIndex(int n, int k);

class VariableGPU: public Variable {
  const int threadsPerBlock_x;
  const int threadsPerBlock_y;

  cufftHandle cufftForwardPlan;
  cufftHandle cufftBackwardPlan;

  public:
    gpu_mode * data_d;
    real * spatialData_d;

    void initialiseData(mode initialValue = 0.0);
    VariableGPU(const Constants &c_in, int totalSteps_in=1, bool useSinTransform_in = true);
    void update(const VariableGPU& dVardt, const real dt, const real f=1.0);
    void fill(const mode value);

    void copyToDevice(bool copySpatial=false);
    void copyToHost(bool copySpatial=false);

    gpu_mode * getCurrent();
    const gpu_mode * getCurrent() const;
    gpu_mode* getPrevious();
    const gpu_mode* getPrevious() const;
    gpu_mode * getPlus(int nSteps=0);
    const gpu_mode* getPlus(int nSteps) const;

    void setupFFTW();
    void toSpectral();
    void toPhysical();
    void postFFTNormalise();

    void readFromFile(std::ifstream& file);
    void writeToFile(std::ofstream& file);
    ~VariableGPU();
};
