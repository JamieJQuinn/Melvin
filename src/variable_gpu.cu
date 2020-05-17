#include <variable_gpu.hpp>
#include <complex_gpu.hpp>
#include <precision.hpp>

#include <iostream>

__global__
void gpu_update(gpu_mode *var, const gpu_mode *dVardt, const gpu_mode *dVardtPrevious, const real dt, const real frac,
    const int nN, const int nZ) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=n_index; n<nN; n+=n_stride) {
    for(int k=k_index; k<nZ; k+=k_stride) {
      int i=k+n*nZ;
      var[i] += ((2.0+frac)*dVardt[i] - frac*dVardtPrevious[i])*dt*0.5;
    }
  }
}

void VariableGPU::initialiseData(real initialValue) {
  cudaMallocManaged(&data, totalSize()*sizeof(real));
  fill(initialValue);
}

void VariableGPU::update(const VariableGPU& dVardt, const real dt, const real f) {
  dim3 threadsPerBlock(threadsPerBlock_x,threadsPerBlock_y);
  dim3 numBlocks((nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_update<<<numBlocks,threadsPerBlock>>>(this->getPlus(), dVardt.getCurrent(), dVardt.getPrevious(), dt, f, nN, nZ);
}

void VariableGPU::readFromFile(std::ifstream& file) {
  real *tempData = new real [totalSize()];
  file.read(reinterpret_cast<char*>(tempData), sizeof(tempData[0])*totalSize());
  for(int i=0; i<totalSteps; ++i) {
    for(int k=0; k<nZ; ++k) {
      for(int n=0; n<nN; ++n) {
        getPlus(i)[k+n*nZ] = makeComplex(tempData[i*nN*nZ + n*nZ + k], 0.0);
      }
    }
  }
  delete [] tempData;
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
  return (gpu_mode*)(data + previous*varSize());
}

const gpu_mode* VariableGPU::getPrevious() const {
  return (gpu_mode*)(data + previous*varSize());
}

gpu_mode* VariableGPU::getPlus(int nSteps) {
  return (gpu_mode*)(data + ((current+nSteps)%totalSteps)*varSize());
}

const gpu_mode* VariableGPU::getPlus(int nSteps) const {
  return (gpu_mode*)(data + ((current+nSteps)%totalSteps)*varSize());
}
