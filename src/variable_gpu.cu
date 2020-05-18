#include <variable_gpu.hpp>
#include <complex_gpu.hpp>
#include <precision.hpp>
#include <gpu_error_checking.hpp>

#include <iostream>

// CUDA constants
__device__ __constant__ int nG_d;
__device__ __constant__ int nX_d;
__device__ __constant__ int nN_d;
__device__ __constant__ int nZ_d;

__device__ int calcIndex(int n, int k) {
  return (k+nG_d)*(nX_d+2*nG_d) + n+nG_d;
}

__global__
void gpu_update(gpu_mode *var, const gpu_mode *dVardt, const gpu_mode *dVardtPrevious, const real dt, const real frac) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=n_index; n<nN_d; n+=n_stride) {
    for(int k=k_index; k<nZ_d; k+=k_stride) {
      int i=calcIndex(n, k);
      var[i] += ((2.0+frac)*dVardt[i] - frac*dVardtPrevious[i])*dt*0.5;
    }
  }
}

void VariableGPU::initialiseData(mode initialValue) {
  Variable::initialiseData(initialValue);
  gpuErrchk(cudaMalloc(&data_d, totalSize()*sizeof(gpu_mode)));
  fill(initialValue);
  gpuErrchk(cudaMemcpyToSymbol(nG_d, &nG, sizeof(nG), 0, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(nX_d, &nX, sizeof(nX), 0, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(nN_d, &nN, sizeof(nN), 0, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(nZ_d, &nZ, sizeof(nZ), 0, cudaMemcpyHostToDevice));
}

void VariableGPU::fill(const mode value) {
  for(int i=0; i<this->totalSize(); ++i) {
    data[i] = value;
  }
  copyToDevice();
}

void VariableGPU::update(const VariableGPU& dVardt, const real dt, const real f) {
  dim3 threadsPerBlock(threadsPerBlock_x,threadsPerBlock_y);
  dim3 numBlocks((nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_update<<<numBlocks,threadsPerBlock>>>(this->getPlus(), dVardt.getCurrent(), dVardt.getPrevious(), dt, f);
}

void VariableGPU::readFromFile(std::ifstream& file) {
  Variable::readFromFile(file);
  copyToDevice();
}

void VariableGPU::writeToFile(std::ofstream& file) {
  copyToHost();
  Variable::writeToFile(file);
}

void VariableGPU::copyToDevice() {
  gpuErrchk(cudaMemcpy(data_d, data, totalSize()*sizeof(data[0]), cudaMemcpyHostToDevice));
}

void VariableGPU::copyToHost() {
  gpuErrchk(cudaMemcpy(data, data_d, totalSize()*sizeof(data[0]), cudaMemcpyDeviceToHost));
}

VariableGPU::VariableGPU(const Constants &c_in, int totalSteps_in, bool useSinTransform_in):
  Variable(c_in, totalSteps_in)
  , threadsPerBlock_x(c_in.threadsPerBlock_x)
  , threadsPerBlock_y(c_in.threadsPerBlock_y)
{}

VariableGPU::~VariableGPU() {
  if(data != NULL) {
    cudaFree(data);
    data = NULL;
  }
}

gpu_mode* VariableGPU::getCurrent() {
  return (gpu_mode*)(getPlus(0));
}

const gpu_mode* VariableGPU::getCurrent() const {
  return (gpu_mode*)(getPlus(0));
}

gpu_mode* VariableGPU::getPrevious() {
  return (gpu_mode*)(data_d + previous*varSize());
}

const gpu_mode* VariableGPU::getPrevious() const {
  return (gpu_mode*)(data_d + previous*varSize());
}

gpu_mode* VariableGPU::getPlus(int nSteps) {
  return (gpu_mode*)(data_d + ((current+nSteps)%totalSteps)*varSize());
}

const gpu_mode* VariableGPU::getPlus(int nSteps) const {
  return (gpu_mode*)(data_d + ((current+nSteps)%totalSteps)*varSize());
}
