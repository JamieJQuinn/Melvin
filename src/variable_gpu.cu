#include <variable_gpu.hpp>
#include <iostream>

__global__
void gpu_update(real *var, const real *dVardt, const real *dVardtPrevious, const real dt, const real frac,
    const int nN, const int nZ) {
  for(int n=0; n<nN; ++n) {
    for(int k=0; k<nZ; ++k) {
      int i=k+n*nZ;
      var[i] += ((2.0+frac)*dVardt[i] - frac*dVardtPrevious[i])*dt*0.5;
    }
  }
}

void VariableGPU::initialiseData(real initialValue) {
  cudaMallocManaged(&data, totalSize()*sizeof(real));
  fill(initialValue);
}

void VariableGPU::update(const Variable& dVardt, const real dt, const real f) {
  gpu_update<<<1,1>>>(this->getPlus(), dVardt.getCurrent(), dVardt.getPrevious(), dt, f,
      nN, nZ);
}

VariableGPU::VariableGPU(const Constants &c_in, int totalSteps_in):
  Variable(c_in, totalSteps_in)
{}

VariableGPU::~VariableGPU() {
  if(data != NULL) {
    cudaFree(data);
    data = NULL;
  }
}
