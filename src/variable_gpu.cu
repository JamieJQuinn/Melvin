#include <variable_gpu.hpp>
#include <complex_gpu.hpp>
#include <precision.hpp>
#include <gpu_error_checking.hpp>

#include <iostream>

// CUDA constants
__device__ __constant__ int extern nG_d;
__device__ __constant__ int extern nX_d;
__device__ __constant__ int extern nN_d;
__device__ __constant__ int extern nZ_d;

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

__global__
void normalise_fft(gpu_mode *data) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=n_index; n<nN_d; n+=n_stride) {
    for(int k=k_index; k<nZ_d; k+=k_stride) {
      int i=calcIndex(n, k);
      data[i] = data[i]*(1.0/nX_d);
    }
  }
}

__global__
void gpu_applyPeriodicHorizontalBoundaryConditions(real *data) {
  int index = threadIdx.x;
  int stride = blockDim.x;
  for(int k=index; k<nZ_d; k+=stride) {
    data[calcIndex(-1, k)] = data[calcIndex(nX_d-1,k)];
    data[calcIndex(nX_d, k)] = data[calcIndex(0,k)];
  }
}

void VariableGPU::applyVerticalBoundaryConditions() {
  if(verticalBoundaryConditions == BoundaryConditions::dirichlet) {
    int stripSize = (nG+1)*rowSize();
    gpuErrchk(cudaMemcpy(data_d, bottomBoundary, stripSize*sizeof(bottomBoundary[0]), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(data_d+varSize()-stripSize, topBoundary, stripSize*sizeof(topBoundary[0]), cudaMemcpyHostToDevice));
  } else if (verticalBoundaryConditions == BoundaryConditions::periodic) {
    // TODO copy top to bottom and bottom to top
  }
}

void VariableGPU::initialiseData(mode initialValue) {
  gpuErrchk(cudaMalloc(&data_d, totalSize()*sizeof(gpu_mode)));
  gpuErrchk(cudaMalloc(&spatialData_d, totalSize()*sizeof(real)));
  fill(initialValue);
}

void VariableGPU::fill(const mode value) {
  for(int i=0; i<this->totalSize(); ++i) {
    data[i] = value;
    spatialData[i] = value.real();
  }
  copyToDevice(true);
}

void VariableGPU::update(const VariableGPU& dVardt, const real dt, const real f) {
  dim3 threadsPerBlock(threadsPerBlock_x,threadsPerBlock_y);
  dim3 numBlocks((nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_update<<<numBlocks,threadsPerBlock>>>(this->getPlus(), dVardt.getCurrent(), dVardt.getPrevious(), dt, f);
}

void VariableGPU::readFromFile(std::ifstream& file) {
  Variable::readFromFile(file);
  copyToDevice();
  fillVerticalBoundaryConditions();
}

void VariableGPU::fillVerticalBoundaryConditions() {
  // Take all ghost rows and one regular row as boundary
  int stripSize = (nG+1)*rowSize();
  for(int i=0; i<stripSize; ++i) {
    bottomBoundary[i] = data[i];
    int dataIdx = varSize()-stripSize+i;
    topBoundary[i] = data[dataIdx];
  }
}

void VariableGPU::applyPhysicalHorizontalBoundaryConditions() {
  if(horizontalBoundaryConditions == BoundaryConditions::periodic) {
    dim3 threadsPerBlock(c.threadsPerBlock_x*c.threadsPerBlock_y);
    dim3 numBlocks((c.nZ - 1 + threadsPerBlock.x)/threadsPerBlock.x);
    gpu_applyPeriodicHorizontalBoundaryConditions<<<numBlocks,threadsPerBlock>>>(spatialData_d);
  }
}

//void VariableGPU::applyPhysicalVerticalBoundaryConditions() {
  //real nX = c.nX;
  //if(c.horizontalBoundaryConditions == BoundaryConditions::impermeable) {
    //for(int k=0; k<c.nZ; ++k) {
      //// Non-conducting
      //vars.tmp.spatial(-1,k) = vars.tmp.spatial(1, k);
      //vars.tmp.spatial(nX,k) = vars.tmp.spatial(nX-2, k);

      //// Impermeable
      //vars.psi.spatial(0,k) = 0.0;
      //vars.psi.spatial(nX-1,k) = 0.0;

      //// Stress free
      //vars.psi.spatial(-1,k) = 2.0*vars.psi.spatial(0,k) - vars.psi.spatial(1,k);
      //vars.psi.spatial(nX,k) = 2.0*vars.psi.spatial(nX-1,k) - vars.psi.spatial(nX-2,k);
      //vars.omg.spatial(0,k) = 0.0;
      //vars.omg.spatial(nX-1,k) = 0.0;
    //}
  //} else if(c.horizontalBoundaryConditions == BoundaryConditions::periodic) {
    //for(int k=0; k<c.nZ; ++k) {
      //vars.tmp.spatial(-1,k) = vars.tmp.spatial(nX-1, k);
      //vars.tmp.spatial(nX,k) = vars.tmp.spatial(0, k);

      //vars.omg.spatial(-1,k) = vars.omg.spatial(nX-1, k);
      //vars.omg.spatial(nX,k) = vars.omg.spatial(0, k);

      //vars.psi.spatial(-1,k) = vars.psi.spatial(nX-1, k);
      //vars.psi.spatial(nX,k) = vars.psi.spatial(0, k);
    //}
  //}
//}

void VariableGPU::writeToFile(std::ofstream& file) {
  copyToHost();
  Variable::writeToFile(file);
}

void VariableGPU::copyToDevice(bool copySpatial) {
  gpuErrchk(cudaMemcpy(data_d, data, totalSize()*sizeof(data[0]), cudaMemcpyHostToDevice));
  if(copySpatial) {
    gpuErrchk(cudaMemcpy(spatialData_d, spatialData, totalSize()*sizeof(spatialData[0]), cudaMemcpyHostToDevice));
  }
}

void VariableGPU::copyToHost(bool copySpatial) {
  gpuErrchk(cudaMemcpy(data, data_d, totalSize()*sizeof(data[0]), cudaMemcpyDeviceToHost));
  if(copySpatial) {
    gpuErrchk(cudaMemcpy(spatialData, spatialData_d, totalSize()*sizeof(spatialData[0]), cudaMemcpyDeviceToHost));
  }
}

VariableGPU::VariableGPU(const Constants &c_in, int totalSteps_in, bool useSinTransform_in):
  Variable(c_in, totalSteps_in)
  , data_d(nullptr)
  , spatialData_d(nullptr)
  , threadsPerBlock_x(c_in.threadsPerBlock_x)
  , threadsPerBlock_y(c_in.threadsPerBlock_y)
{
  initialiseData();
  setupcuFFT();
  topBoundary = new mode[(nG+1)*rowSize()];
  bottomBoundary = new mode[(nG+1)*rowSize()];
}

VariableGPU::~VariableGPU() {
  if(data_d != nullptr) {
    cudaFree(data_d);
    data = nullptr;
  }
  if(spatialData_d != nullptr) {
    cudaFree(spatialData_d);
    data = nullptr;
  }
  delete [] topBoundary;
  delete [] bottomBoundary;
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

void VariableGPU::setupcuFFT() {
  int rank = 1;
  int n[] = {nX};
  int inembed[] = {rowSize()};
  int istride = 1;
  int idist = rowSize();
  int onembed[] = {rowSize()};
  int ostride = 1;
  int odist = rowSize();
  int batch = nZ;

  cufftType type = CUFFT_D2Z;
  cufftResult result = cufftPlanMany(&cufftForwardPlan,
      rank, n,
      inembed, istride, idist,
      onembed, ostride, odist,
      type, batch);

  if(result != CUFFT_SUCCESS) {
    std::cerr << "cuFFT forward plan could not be created" << std::endl;
  }

  type = CUFFT_Z2D;
  result = cufftPlanMany(&cufftBackwardPlan,
      rank, n,
      inembed, istride, idist,
      onembed, ostride, odist,
      type, batch);

  if(result != CUFFT_SUCCESS) {
    std::cerr << "cuFFT backward plan could not be created" << std::endl;
  }
}

void VariableGPU::postFFTNormalise() {
  dim3 threadsPerBlock(threadsPerBlock_x,threadsPerBlock_y);
  dim3 numBlocks((nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (nZ + threadsPerBlock.y - 1)/threadsPerBlock.y);
  normalise_fft<<<numBlocks,threadsPerBlock>>>(getCurrent());
}

void VariableGPU::toSpectral() {
  real *spatial = spatialData_d + calcIndex(0,0);
  gpu_mode *spectral = getCurrent() + calcIndex(0,0);

  cufftResult result = cufftExecD2Z(cufftForwardPlan, (cufftDoubleReal*)spatial, (cufftDoubleComplex*)spectral);

  if(result != CUFFT_SUCCESS) {
    std::cerr << "cuFFT forward plan could not be executed. Error:" << result << std::endl;
  }

  postFFTNormalise();
  cudaDeviceSynchronize();
}

void VariableGPU::toPhysical() {
  real *spatial = spatialData_d + calcIndex(0,0);
  gpu_mode *spectral = getCurrent() + calcIndex(0,0);

  cufftResult result = cufftExecZ2D(cufftBackwardPlan, (cufftDoubleComplex*)spectral, (cufftDoubleReal*)spatial);
  if(result != CUFFT_SUCCESS) {
    std::cerr << "cuFFT backward plan could not be executed. Error:"<< result << std::endl;
  }
}
