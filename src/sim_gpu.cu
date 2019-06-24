#include <sim_gpu.hpp>

#include <precision.hpp>
#include <math.h>

__device__
real dfdz2(const real *data, const int n, const int k, const int nZ, const int oodz2) {
  int in = n*nZ;

  return (data[k+1 + in] - 2.0f*data[k + in] + data[k-1 + in])*oodz2;
}

__device__
real dfdz(const real *data, const int n, const int k, const int nZ, const int oodz) {
  int in = n*nZ;

  return (data[k+1 + in] - data[k-1 + in])*oodz*0.5f;
}

__device__
real sqr(real x) {
  return x*x;
}

__global__
void gpu_computeLinearTemperatureDerivative(real *dTmpdt, const real *tmp,
    const int nN, const int nZ, const real aspectRatio, const real oodz2) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=n_index; n<nN; n+=n_stride) {
    for(int k=1+k_index; k<nZ-1; k+=k_stride) {
      int i=k+n*nZ;
      dTmpdt[i] = dfdz2(tmp, n, k, nZ, oodz2) - sqr(n*M_PI/aspectRatio)*tmp[i];
    }
  }
}

__global__
void gpu_computeLinearVorticityDerivative(real *dOmgdt, const real *omg, const real *tmp,
    const int nN, const int nZ, const real aspectRatio, const real Ra, const real Pr, const real oodz2) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=n_index; n<nN; n+=n_stride) {
    for(int k=1+k_index; k<nZ-1; k+=k_stride) {
      int i=k+n*nZ;
      dOmgdt[i] =
        Pr*(dfdz2(omg,n,k,nZ,oodz2) - sqr(n*M_PI/aspectRatio)*omg[i])
        + Pr*Ra*(n*M_PI/aspectRatio)*tmp[i];
    }
  }
}

__global__
void gpu_computeLinearXiDerivative(real *dXidt, const real *xi, real *dOmgdt, const real *omg,
    const int nN, const int nZ, const real tau, const real aspectRatio, const real RaXi, const real Pr, const real oodz2) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=n_index; n<nN; n+=n_stride) {
    for(int k=1+k_index; k<nZ-1; k+=k_stride) {
      int i=k+n*nZ;
      dXidt[i] = tau*(dfdz2(xi, n, k, nZ, oodz2) - pow(n*M_PI/aspectRatio, 2)*xi[i]);
      dOmgdt[i] += -RaXi*tau*Pr*(n*M_PI/aspectRatio)*xi[i];
    }
  }
}

__global__
void gpu_fillMode(real *data, const real value, const int n, const int nZ) {
  int index = threadIdx.x;
  int stride = blockDim.x;
  for(int k=index; k<nZ; k+=stride) {
    int i=k+n*nZ;
    data[i] = value;
  }
}

__global__
void gpu_addAdvectionApproximation(
    real *dVardt, const real *var,
    const real *psi,
    const int nN, const int nZ, const real aspectRatio, const real oodz) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=1+n_index; n<nN; n+=n_stride) {
    for(int k=1+k_index; k<nZ-1; k+=k_stride) {
      int i=k+n*nZ;
      dVardt[i] += -1*dfdz(var,0,k,nZ,oodz)*n*M_PI/aspectRatio * psi[i];
    }
  }
}

SimGPU::SimGPU(const Constants &c_in)
  : c(c_in)
  , vars(c_in)
  , keTracker(c_in)
{
  dt = c.initialDt;

  thomasAlgorithm = new ThomasAlgorithmGPU(c.nZ, c.nN, c.aspectRatio, c.oodz2);
}

SimGPU::~SimGPU() {
  delete thomasAlgorithm;
}

void SimGPU::computeLinearTemperatureDerivative() {
  dim3 threadsPerBlock(8,16);
  dim3 numBlocks((c.nN + 8 - 1)/8, (c.nZ - 2 + 16 - 1)/16);
  gpu_computeLinearTemperatureDerivative<<<numBlocks,threadsPerBlock>>>(vars.dTmpdt.getCurrent(), vars.tmp.getCurrent(), c.nN, c.nZ, c.aspectRatio, c.oodz2);
}

void SimGPU::computeLinearVorticityDerivative() {
  dim3 threadsPerBlock(8,16);
  dim3 numBlocks((c.nN + 8 - 1)/8, (c.nZ - 2 + 16 - 1)/16);
  gpu_computeLinearVorticityDerivative<<<numBlocks,threadsPerBlock>>>(vars.dOmgdt.getCurrent(), vars.omg.getCurrent(), vars.tmp.getCurrent(),
    c.nN, c.nZ, c.aspectRatio, c.Ra, c.Pr, c.oodz2);
}

void SimGPU::computeLinearXiDerivative() {
  dim3 threadsPerBlock(8,16);
  dim3 numBlocks((c.nN + 8 - 1)/8, (c.nZ - 2 + 16 - 1)/16);
  gpu_computeLinearXiDerivative<<<numBlocks,threadsPerBlock>>>(vars.dXidt.getCurrent(), vars.xi.getCurrent(), vars.dOmgdt.getCurrent(), vars.omg.getCurrent(),
    c.nN, c.nZ, c.tau, c.aspectRatio, c.RaXi, c.Pr, c.oodz2);
}

void SimGPU::computeLinearDerivatives() {
  // Computes the (linear) derivatives of Tmp and vars.omg
  computeLinearTemperatureDerivative();
  computeLinearVorticityDerivative();
  if(c.isDoubleDiffusion) {
    computeLinearXiDerivative();
  }
}

void SimGPU::addAdvectionApproximation() {
  // Only applies to the linear simulation
  dim3 threadsPerBlock(8,16);
  dim3 numBlocks((c.nN -1 + 8 - 1)/8, (c.nZ + 16 - 1)/16);
  gpu_fillMode<<<1,256>>>(vars.dOmgdt.getCurrent(), 0.0, 0, c.nZ);
  gpu_fillMode<<<1,256>>>(vars.dTmpdt.getCurrent(), 0.0, 0, c.nZ);
  gpu_addAdvectionApproximation<<<numBlocks,threadsPerBlock>>>(
      vars.dTmpdt.getCurrent(), vars.tmp.getCurrent(), vars.psi.getCurrent(),
      c.nN, c.nZ, c.aspectRatio, 1.0f/c.dz);
  if(c.isDoubleDiffusion) {
    gpu_fillMode<<<1,256>>>(vars.dXidt.getCurrent(), 0.0, 0, c.nZ);
    gpu_addAdvectionApproximation<<<1,1>>>(
        vars.dXidt.getCurrent(), vars.xi.getCurrent(), vars.psi.getCurrent(),
        c.nN, c.nZ, c.aspectRatio, 1.0f/c.dz);
  }
}

void SimGPU::solveForPsi(){
  // Solve for Psi using Thomas algorithm
  thomasAlgorithm->solve(vars.psi.getCurrent(), vars.omg.getCurrent());
}

void SimGPU::runLinearStep() {
  computeLinearDerivatives();
  addAdvectionApproximation();
  vars.updateVars(dt);
  vars.advanceDerivatives();
  solveForPsi();
}
