#include <sim_gpu.hpp>
#include <variable_gpu.hpp>

#include <precision.hpp>
#include <numerical_methods.hpp>
#include <complex_gpu.hpp>
#include <gpu_error_checking.hpp>

#include <math.h>
#include <iostream>

using std::cout;
using std::endl;

__device__ __constant__ extern int nX_d;
__device__ __constant__ extern int nN_d;
__device__ __constant__ extern int nZ_d;
__device__ __constant__ extern real oodz_d;
__device__ __constant__ extern real oodx_d;
__device__ __constant__ extern real oodz2_d;
__device__ __constant__ extern real aspectRatio_d;
__device__ __constant__ extern real wavelength_d;
__device__ __constant__ extern gpu_mode xSinDerivativeFactor_d;
__device__ __constant__ extern gpu_mode xCosDerivativeFactor_d;
__device__ __constant__ extern real Ra_d;
__device__ __constant__ extern real Pr_d;
__device__ __constant__ extern real RaXi_d;
__device__ __constant__ extern real tau_d;

__device__
gpu_mode dfdz2(const gpu_mode *data, const int n, const int k) {
  return (data[calcIndex(n, k+1)] - 2.0f*data[calcIndex(n, k)] + data[calcIndex(n, k-1)])*oodz2_d;
}

__device__
gpu_mode dfdz(const gpu_mode *data, const int n, const int k) {
  return (data[calcIndex(n, k+1)] - data[calcIndex(n, k-1)])*oodz_d*0.5;
}

__device__
real dfdz_s(const real *data, const int i, const int k) {
  return (data[calcIndex(i, k+1)] - data[calcIndex(i, k-1)])*oodz_d*0.5;
}

__device__
real dfdx_s(const real *data, const int i, const int k) {
  return (data[calcIndex(i+1, k)] - data[calcIndex(i-1, k)])*oodx_d*0.5;
}

__device__
gpu_mode sqr(gpu_mode x) {
  return x*x;
}

__device__
gpu_mode sqr(real x) {
  return sqr(makeComplex(x, 0.0));
}

__global__
void gpu_computeLinearTemperatureDerivative(gpu_mode *dTmpdt, const gpu_mode *tmp) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=n_index; n<nN_d; n+=n_stride) {
    for(int k=k_index; k<nZ_d; k+=k_stride) {
      int i=calcIndex(n,k);
      dTmpdt[i] = dfdz2(tmp, n, k) - sqr(n*wavelength_d)*tmp[i];
    }
  }
}

__global__
void gpu_computeLinearVorticityDerivative(gpu_mode *dOmgdt, const gpu_mode *omg, const gpu_mode *tmp) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=n_index; n<nN_d; n+=n_stride) {
    for(int k=k_index; k<nZ_d; k+=k_stride) {
      int i=calcIndex(n,k);
      dOmgdt[i] =
        Pr_d*(dfdz2(omg,n,k) - sqr(n*wavelength_d)*omg[i])
        - xCosDerivativeFactor_d*Pr_d*Ra_d*(n*wavelength_d)*tmp[i];
    }
  }
}

__global__
void gpu_computeLinearXiDerivative(gpu_mode *dXidt, const gpu_mode *xi, gpu_mode *dOmgdt, const gpu_mode *omg) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=n_index; n<nN_d; n+=n_stride) {
    for(int k=k_index; k<nZ_d; k+=k_stride) {
      int i=calcIndex(n,k);
      dXidt[i] = tau_d*(dfdz2(xi, n, k) - sqr(n*wavelength_d)*xi[i]);
      dOmgdt[i] += xCosDerivativeFactor_d*RaXi_d*tau_d*Pr_d*(n*wavelength_d)*xi[i];
    }
  }
}

__global__
void gpu_fillMode(gpu_mode *data, const gpu_mode value, const int n) {
  int index = threadIdx.x;
  int stride = blockDim.x;
  for(int k=index; k<nZ_d; k+=stride) {
    int i=calcIndex(n,k);
    data[i] = value;
  }
}

__global__
void gpu_addAdvectionApproximation(
    gpu_mode *dVardt, const gpu_mode *var,
    const gpu_mode *psi) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=1+n_index; n<nN_d; n+=n_stride) {
    for(int k=k_index; k<nZ_d; k+=k_stride) {
      int i=calcIndex(n,k);
      dVardt[i] += -1*xSinDerivativeFactor_d*dfdz(var,0,k)*(n*wavelength_d)*psi[i];
    }
  }
}

__global__
void gpu_plusEquals(
    gpu_mode *var1, const gpu_mode *var2) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=n_index; n<nN_d; n+=n_stride) {
    for(int k=k_index; k<nZ_d; k+=k_stride) {
      int i=calcIndex(n,k);
      var1[i] += var2[i];
    }
  }
}


__global__
void gpu_computeNonlinearDerivativeN0(
    gpu_mode *dVardt, const gpu_mode *var,
    const gpu_mode *psi) {
  int k_index = blockIdx.x*blockDim.x + threadIdx.x;
  int k_stride = blockDim.x*gridDim.x;
  for(int k=k_index; k<nZ_d; k+=k_stride) {
    for(int n=1; n<nN_d; ++n) {
      // Contribution TO var[n=0]
      int i=calcIndex(n,k);
      dVardt[calcIndex(0,k)] +=
        -M_PI/(2*aspectRatio_d)*n*(
          dfdz(psi,n,k)*var[i] +
          dfdz(var,n,k)*psi[i]
          );
      // Contribution FROM var[n=0]
      dVardt[i] +=
        -n*M_PI/aspectRatio_d*psi[i]*dfdz(var,0,k);
    }
  }
}

__global__
void gpu_computeNonlinearDerivative(
    gpu_mode *dVardt, const gpu_mode *var,
    const gpu_mode *psi,
    const int vorticityFactor) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int k=k_index; k<nZ_d; k+=k_stride) {
    for(int n=1+n_index; n<nN_d; n+=n_stride) {
      // Contribution FROM var[n>0] and vars.omg[n>0]
      int o;
      for(int m=1; m<n; ++m){
        // Case n = n' + n''
        o = n-m;
        int im = calcIndex(m,k);
        int io = calcIndex(o,k);
        int in = calcIndex(n,k);
        dVardt[in] +=
          -M_PI/(2.0*aspectRatio_d)*(
          -m*dfdz(psi,o,k)*var[im]
          +o*dfdz(var,m,k)*psi[io]
          );
      }
      for(int m=n+1; m<nN_d; ++m){
        // Case n = n' - n''
        o = m-n;
        int im = calcIndex(m,k);
        int io = calcIndex(o,k);
        int in = calcIndex(n,k);
        dVardt[in] +=
          -M_PI/(2.0*aspectRatio_d)*(
          +m*dfdz(psi,o,k)*var[im]
          +o*dfdz(var,m,k)*psi[io]
          );
      }
      for(int m=1; m+n<nN_d; ++m){
        // Case n= n'' - n'
        o = n+m;
        int im = calcIndex(m,k);
        int io = calcIndex(o,k);
        int in = calcIndex(n,k);
        dVardt[in] +=
          vorticityFactor*M_PI/(2.0*aspectRatio_d)*(
          +m*dfdz(psi,o,k)*var[im]
          +o*dfdz(var,m,k)*psi[io]
          );
      }
    }
  }
}

__global__
void gpu_computeNonlinearDerivativeSpectralTransform(real *nonlinearTerm, const real *var, const real *psi) {
  int i_index = blockIdx.x*blockDim.x + threadIdx.x;
  int i_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int i=i_index; i<nX_d; i+=i_stride) {
    for(int k=k_index; k<nZ_d; k+=k_stride) {
      int ix = calcIndex(i, k);
      nonlinearTerm[ix] = 
        -(
            (var[calcIndex(i+1,k)]*(-dfdz_s(psi,i+1,k)) -
             var[calcIndex(i-1,k)]*(-dfdz_s(psi,i-1,k)))*oodx_d*0.5 +
            (var[calcIndex(i,k+1)]* dfdx_s(psi,i,k+1) -
             var[calcIndex(i,k-1)]* dfdx_s(psi,i,k-1))*oodz_d*0.5
         );
    }
  }
}

SimGPU::SimGPU(const Constants &c_in)
  : c(c_in)
  , vars(c_in)
  , keTracker(c_in)
  , nonlinearTerm(c_in)
  , threadsPerBlock2D(c_in.threadsPerBlock_x, c_in.threadsPerBlock_y)
  , numBlocks2D(
      (c_in.nN + c_in.threadsPerBlock_x - 1)/c_in.threadsPerBlock_x,
      (c_in.nZ + c_in.threadsPerBlock_y - 1)/c_in.threadsPerBlock_y
    )
  , numBlocks2DSpatial(
      (c_in.nX + c_in.threadsPerBlock_x - 1)/c_in.threadsPerBlock_x,
      (c_in.nZ + c_in.threadsPerBlock_y - 1)/c_in.threadsPerBlock_y
    )
  , threadsPerBlock1D(c_in.threadsPerBlock_x*c_in.threadsPerBlock_y)
  , numBlocks1D(
      (c_in.nZ + c_in.threadsPerBlock_y - 1)/c_in.threadsPerBlock_y
    )
{
  dt = c.initialDt;

  thomasAlgorithm = new ThomasAlgorithmGPU(c);
}

SimGPU::~SimGPU() {
  delete thomasAlgorithm;
}

void SimGPU::computeLinearTemperatureDerivative() {
  gpu_computeLinearTemperatureDerivative<<<numBlocks2D,threadsPerBlock2D>>>(vars.dTmpdt.getCurrent(), vars.tmp.getCurrent());
}

void SimGPU::computeLinearVorticityDerivative() {
  gpu_computeLinearVorticityDerivative<<<numBlocks2D,threadsPerBlock2D>>>(vars.dOmgdt.getCurrent(), vars.omg.getCurrent(), vars.tmp.getCurrent());
}

void SimGPU::computeLinearXiDerivative() {
  gpu_computeLinearXiDerivative<<<numBlocks2D,threadsPerBlock2D>>>(vars.dXidt.getCurrent(), vars.xi.getCurrent(), vars.dOmgdt.getCurrent(), vars.omg.getCurrent());
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
  gpu_fillMode<<<numBlocks1D,threadsPerBlock1D>>>(vars.dOmgdt.getCurrent(), makeComplex(0.0, 0.0), 0);
  gpu_fillMode<<<numBlocks1D,threadsPerBlock1D>>>(vars.dTmpdt.getCurrent(), makeComplex(0.0, 0.0), 0);
  gpu_addAdvectionApproximation<<<numBlocks2D,threadsPerBlock2D>>>(
      vars.dTmpdt.getCurrent(), vars.tmp.getCurrent(), vars.psi.getCurrent());
  if(c.isDoubleDiffusion) {
    gpu_fillMode<<<numBlocks1D,threadsPerBlock1D>>>(vars.dXidt.getCurrent(), makeComplex(0.0, 0.0), 0);
    gpu_addAdvectionApproximation<<<numBlocks2D,threadsPerBlock2D>>>(
        vars.dXidt.getCurrent(), vars.xi.getCurrent(), vars.psi.getCurrent());
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
  vars.tmp.applyVerticalBoundaryConditions();
  vars.omg.applyVerticalBoundaryConditions();
  vars.advanceDerivatives();
  solveForPsi();
  vars.psi.applyVerticalBoundaryConditions();
}

void SimGPU::computeNonlinearDerivativeSpectralTransform(VariableGPU& dVardt, const VariableGPU& var) {
  gpu_computeNonlinearDerivativeSpectralTransform<<<numBlocks2DSpatial,threadsPerBlock2D>>>(nonlinearTerm.spatialData_d, var.spatialData_d, vars.psi.spatialData_d);

  nonlinearTerm.toSpectral();

  gpu_plusEquals<<<numBlocks2D,threadsPerBlock2D>>>(dVardt.getCurrent(), nonlinearTerm.getCurrent());
}

void SimGPU::computeNonlinearTemperatureDerivative() {
  // Calculate n=0 gpu_mode
  gpu_computeNonlinearDerivativeN0<<<numBlocks1D,threadsPerBlock1D>>>(vars.dTmpdt.getCurrent(), vars.tmp.getCurrent(), vars.psi.getCurrent());

  // Calculate other gpu_modes
  gpu_computeNonlinearDerivative<<<numBlocks2D,threadsPerBlock2D>>>(vars.dTmpdt.getCurrent(), vars.tmp.getCurrent(), vars.psi.getCurrent(), -1);
}

void SimGPU::computeNonlinearXiDerivative() {
  // Calculate n=0 gpu_mode
  gpu_computeNonlinearDerivativeN0<<<numBlocks1D,threadsPerBlock1D>>>(vars.dXidt.getCurrent(), vars.xi.getCurrent(), vars.psi.getCurrent());

  // Calculate other gpu_modes
  gpu_computeNonlinearDerivative<<<numBlocks2D,threadsPerBlock2D>>>(vars.dXidt.getCurrent(), vars.xi.getCurrent(), vars.psi.getCurrent(), -1);
}

void SimGPU::computeNonlinearVorticityDerivative() {
  gpu_computeNonlinearDerivative<<<numBlocks2D,threadsPerBlock2D>>>(vars.dOmgdt.getCurrent(), vars.omg.getCurrent(), vars.psi.getCurrent(), 1);
}

void SimGPU::runNonLinearStep(real f) {
  computeLinearDerivatives();
  computeNonlinearDerivatives();
  vars.updateVars(dt, f);
  vars.tmp.applyVerticalBoundaryConditions();
  vars.omg.applyVerticalBoundaryConditions();
  if(c.isDoubleDiffusion) {
    vars.xi.applyVerticalBoundaryConditions();
  }
  vars.advanceDerivatives();
  solveForPsi();
  vars.psi.applyVerticalBoundaryConditions();
}

void SimGPU::computeNonlinearDerivativesGalerkin() {
  computeNonlinearTemperatureDerivative();
  computeNonlinearVorticityDerivative();
  if(c.isDoubleDiffusion) {
    computeNonlinearXiDerivative();
  }
}

void SimGPU::computeNonlinearDerivativesSpectralTransform() {
  vars.psi.toPhysical();
  vars.tmp.toPhysical();
  vars.omg.toPhysical();

  vars.psi.applyPhysicalHorizontalBoundaryConditions();
  vars.tmp.applyPhysicalHorizontalBoundaryConditions();
  vars.omg.applyPhysicalHorizontalBoundaryConditions();

  computeNonlinearDerivativeSpectralTransform(vars.dTmpdt, vars.tmp);
  computeNonlinearDerivativeSpectralTransform(vars.dOmgdt, vars.omg);
  if(c.isDoubleDiffusion) {
    vars.xi.toPhysical();
    vars.xi.applyPhysicalHorizontalBoundaryConditions();
    computeNonlinearDerivativeSpectralTransform(vars.dXidt, vars.xi);
  }
}

void SimGPU::computeNonlinearDerivatives() {
  if(c.horizontalBoundaryConditions == BoundaryConditions::periodic) {
    computeNonlinearDerivativesSpectralTransform();
  } else {
    computeNonlinearDerivativesGalerkin();
  }
}

void SimGPU::runNonLinear() {
  vars.load(c.icFile);

  vars.psi.applyVerticalBoundaryConditions();
  vars.tmp.applyVerticalBoundaryConditions();
  vars.omg.applyVerticalBoundaryConditions();
  if(c.isDoubleDiffusion) {
    vars.xi.applyVerticalBoundaryConditions();
  }

  real saveTime = 0;
  real KEcalcTime = 0;
  real KEsaveTime = 0;
  real CFLCheckTime = 0;
  real f = 1.0f; // Fractional change in dt (if CFL condition being breached)
  t = 0;
  while (c.totalTime-t>EPSILON) {
    //if(KEcalcTime-t < EPSILON) {
      //cudaDeviceSynchronize();
      //keTracker.calcKineticEnergy(vars.psi);
      //KEcalcTime += 1e2*dt;
    //}
    //if(KEsaveTime-t < EPSILON) {
      //keTracker.saveKineticEnergy();
      //KEsaveTime += 1e4*dt;
    //}
    //if(CFLCheckTime-t < EPSILON) {
      //cout << "Checking CFL" << endl;
      //CFLCheckTime += 1e4*dt;
      //cudaDeviceSynchronize();
      //f = checkCFL(vars.psi, c.dz, c.dx, dt, c.aspectRatio, c.nN, c.nX, c.nZ);
      //dt*=f;
    //}
    if(saveTime-t < EPSILON) {
      cout << t << " of " << c.totalTime << "(" << t/c.totalTime*100 << "%)" << endl;
      saveTime+=c.timeBetweenSaves;
      cudaDeviceSynchronize();
      vars.save();
    }
    runNonLinearStep(f);
    t+=dt;
    f=1.0f;
  }
  cudaDeviceSynchronize();
  printf("%e of %e (%.2f%%)\n", t, c.totalTime, t/c.totalTime*100);
  vars.save();
  //keTracker.saveKineticEnergy();
}
