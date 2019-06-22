#include <variable_gpu.hpp>
#include <iostream>

void VariableGPU::initialiseData(real initialValue) {
  cudaMallocManaged(&data, totalSize()*sizeof(real));
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
