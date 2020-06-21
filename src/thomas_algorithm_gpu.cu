#include "thomas_algorithm_gpu.hpp"
#include <variable_gpu.hpp>
#include <gpu_error_checking.hpp>

__global__
void solveThomasAlgorithm(gpu_mode *sol, const gpu_mode *rhs, const real *wk1, const real *wk2, const real *sub, const int nN, const int nZ) {
  int mode = threadIdx.x;
  int stride = blockDim.x;
  for(int n=mode; n<nN; n+=stride) {
    int iN = n*nZ;

    // Forward Subsitution
    sol[calcIndex(n,0)] = rhs[calcIndex(n,0)]*wk1[0+iN];
    for (int k=1; k<nZ; ++k) {
      sol[calcIndex(n,k)] = (rhs[calcIndex(n,k)] - sub[k-1]*sol[calcIndex(n,k-1)])*wk1[k+iN];
    }
    // Backward Substitution
    for (int k=nZ-2; k>=0; --k) {
      sol[calcIndex(n,k)] -= wk2[k+iN]*sol[calcIndex(n,k+1)];
    }
  }
}

void ThomasAlgorithmGPU::solve(gpu_mode *sol, const gpu_mode *rhs) const {
  solveThomasAlgorithm<<<1,256>>>((gpu_mode*)sol, (gpu_mode*)rhs, (real*)wk1, (real*)wk2, (real*)sub, nN, nZ);
}

void ThomasAlgorithmGPU::formTriDiagonalArraysForN (
          const real *sub, const real *dia, const real *sup,
    real * wk1, real *wk2) {

  wk1[0] = 1.0/dia[0];
  wk2[0] = sup[0]*wk1[0];

  for (int i=1; i<nZ-1; ++i) {
    wk1[i] = 1.0/(dia[i] - sub[i-1] * wk2[i-1]);
    wk2[i] = sup[i]*wk1[i];
  }

  wk1[nZ-1] = 1.0/(dia[nZ-1] - sub[nZ-2]*wk2[nZ-2]);
}

ThomasAlgorithmGPU::ThomasAlgorithmGPU(const Constants& c):
  nZ(c.nZ),
  nN(c.nN),
  oodz2(c.oodz2),
  wavelength(c.wavelength)
  {
  // Precalculate tridiagonal arrays
  real * dia = new real [nZ];
  real * sup = new real [nZ];

  real * subHost = new real [nZ];
  real * wk1Host = new real [nZ*nN];
  real * wk2Host = new real [nZ*nN];

  for(int k=0; k<nZ; ++k) {
    subHost[k] = sup[k] = -oodz2;
  }
  for(int n=0; n<nN; ++n) {
    for(int k=0; k<nZ; ++k){
      dia[k] = pow(wavelength*n, 2) + 2.0*oodz2;
    }
    dia[0] = dia[nZ-1] = 1.0;
    subHost[nZ-2] = sup[0] = 0.0;
    formTriDiagonalArraysForN(
    subHost, dia, sup,
    wk1Host+n*nZ, wk2Host+n*nZ);
  }

  gpuErrchk(cudaMalloc(&sub, nZ*sizeof(real)));
  gpuErrchk(cudaMalloc(&wk1, nZ*nN*sizeof(real)));
  gpuErrchk(cudaMalloc(&wk2, nZ*nN*sizeof(real)));

  gpuErrchk(cudaMemcpy(sub, subHost, nZ*sizeof(real), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(wk1, wk1Host, nZ*nN*sizeof(real), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(wk2, wk2Host, nZ*nN*sizeof(real), cudaMemcpyHostToDevice));

  delete [] dia;
  delete [] sup;

  delete [] subHost;
  delete [] wk1Host;
  delete [] wk2Host;
}

ThomasAlgorithmGPU::~ThomasAlgorithmGPU() {
  cudaFree(wk1);
  cudaFree(wk2);
  cudaFree(sub);
}
