#include <constants_gpu.hpp>
#include <gpu_error_checking.hpp>

__device__ __constant__ int nG_d;
__device__ __constant__ int nX_d;
__device__ __constant__ int nN_d;
__device__ __constant__ int nZ_d;
__device__ __constant__ real oodz_d;
__device__ __constant__ real oodx_d;
__device__ __constant__ real oodz2_d;
__device__ __constant__ real aspectRatio_d;
__device__ __constant__ real Ra_d;
__device__ __constant__ real Pr_d;
__device__ __constant__ real RaXi_d;
__device__ __constant__ real tau_d;

void copyConstantToGPU(const int &hostConstant, int &deviceConstant) {
  gpuErrchk(cudaMemcpyToSymbol(deviceConstant, &hostConstant, sizeof(hostConstant), 0, cudaMemcpyHostToDevice));
}

void copyConstantToGPU(const real &hostConstant, real &deviceConstant) {
  gpuErrchk(cudaMemcpyToSymbol(deviceConstant, &hostConstant, sizeof(hostConstant), 0, cudaMemcpyHostToDevice));
}


void copyGPUConstants(
    int nG, int nX, int nN, int nZ,
    real oodz, real oodx, real oodz2,
    real aspectRatio,
    real Ra, real Pr, real RaXi, real tau
  ) {
  copyConstantToGPU(nG, nG_d);
  copyConstantToGPU(nX, nX_d);
  copyConstantToGPU(nN, nN_d);
  copyConstantToGPU(nZ, nZ_d);
  copyConstantToGPU(oodz, oodz_d);
  copyConstantToGPU(oodx, oodx_d);
  copyConstantToGPU(oodz2, oodz2_d);
  copyConstantToGPU(aspectRatio, aspectRatio_d);
  copyConstantToGPU(Ra, Ra_d);
  copyConstantToGPU(Pr, Pr_d);
  copyConstantToGPU(tau, tau_d);
  copyConstantToGPU(RaXi, RaXi_d);
}
