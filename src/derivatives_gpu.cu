#include <derivatives_gpu.hpp>

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
  for(int n=0; n<nN; ++n) {
    for(int k=1; k<nZ-1; ++k) {
      int i=k+n*nZ;
      dTmpdt[i] = dfdz2(tmp, n, k, nZ, oodz2) - sqr(n*M_PI/aspectRatio)*tmp[i];
    }
  }
}

__global__
void gpu_computeLinearVorticityDerivative(real *dOmgdt, const real *omg, const real *tmp,
    const int nN, const int nZ, const real aspectRatio, const real Ra, const real Pr, const real oodz2) {
  for(int n=0; n<nN; ++n) {
    for(int k=1; k<nZ-1; ++k) {
      int i=k+n*nZ;
      dOmgdt[i] =
        Pr*(dfdz2(omg,n,k,nZ,oodz2) - sqr(n*M_PI/aspectRatio)*omg[i])
        + Pr*Ra*(n*M_PI/aspectRatio)*tmp[i];
    }
  }
}

__global__
void gpu_fillMode(real *data, const real value, const int n, const int nZ) {
  for(int k=0; k<nZ; ++k) {
    int i=k+n*nZ;
    data[i] = value;
  }
}

__global__
void gpu_addAdvectionApproximation(
    real *dVardt, const real *var,
    const real *psi,
    const int nN, const int nZ, const real aspectRatio, const real oodz) {
  for(int n=1; n<nN; ++n) {
    for(int k=1; k<nZ-1; ++k) {
      int i=k+n*nZ;
      dVardt[i] += -1*dfdz(var,0,k,nZ,oodz)*n*M_PI/aspectRatio * psi[i];
    }
  }
}

void computeLinearTemperatureDerivativeGPU(VariableGPU &dTmpdt, const VariableGPU &tmp, const Constants &c) {
  gpu_computeLinearTemperatureDerivative<<<1,1>>>(dTmpdt.getCurrent(), tmp.getCurrent(), c.nN, c.nZ, c.aspectRatio, c.oodz2);
}

void computeLinearVorticityDerivativeGPU(VariableGPU &dOmgdt, const VariableGPU &omg, const VariableGPU &tmp, const Constants &c) {
  gpu_computeLinearVorticityDerivative<<<1,1>>>(dOmgdt.getCurrent(), omg.getCurrent(), tmp.getCurrent(),
    c.nN, c.nZ, c.aspectRatio, c.Ra, c.Pr, c.oodz2);
}

void addAdvectionApproximationGPU(
    VariableGPU &dTmpdt, const VariableGPU &tmp,
    VariableGPU &dOmgdt, const VariableGPU &omg,
    VariableGPU &dXidt, const VariableGPU &xi,
    const VariableGPU &psi,
    const Constants &c) {
  // Only applies to the linear simulation
  gpu_fillMode<<<1,1>>>(dOmgdt.getCurrent(), 0.0, 0, c.nZ);
  gpu_fillMode<<<1,1>>>(dTmpdt.getCurrent(), 0.0, 0, c.nZ);
  gpu_addAdvectionApproximation<<<1,1>>>(
      dTmpdt.getCurrent(), tmp.getCurrent(), psi.getCurrent(),
      c.nN, c.nZ, c.aspectRatio, 1.0f/c.dz);
  if(c.isDoubleDiffusion) {
    gpu_fillMode<<<1,1>>>(dXidt.getCurrent(), 0.0, 0, c.nZ);
    gpu_addAdvectionApproximation<<<1,1>>>(
        dXidt.getCurrent(), xi.getCurrent(), psi.getCurrent(),
        c.nN, c.nZ, c.aspectRatio, 1.0f/c.dz);
  }
}
