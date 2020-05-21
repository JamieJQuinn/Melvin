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
  return (data[calcIndex(i+1, k)] - data[calcIndex(i-1, k)])*oodz_d*0.5;
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
    for(int k=1+k_index; k<nZ_d-1; k+=k_stride) {
      int i=calcIndex(n,k);
      dTmpdt[i] = dfdz2(tmp, n, k) - sqr(n*M_PI/aspectRatio_d)*tmp[i];
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
    for(int k=1+k_index; k<nZ_d-1; k+=k_stride) {
      int i=calcIndex(n,k);
      dOmgdt[i] =
        Pr_d*(dfdz2(omg,n,k) - sqr(n*M_PI/aspectRatio_d)*omg[i])
        + Pr_d*Ra_d*(n*M_PI/aspectRatio_d)*tmp[i];
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
    for(int k=1+k_index; k<nZ_d-1; k+=k_stride) {
      int i=calcIndex(n,k);
      dXidt[i] = tau_d*(dfdz2(xi, n, k) - pow(n*M_PI/aspectRatio_d, 2)*xi[i]);
      dOmgdt[i] += -RaXi_d*tau_d*Pr_d*(n*M_PI/aspectRatio_d)*xi[i];
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
    for(int k=1+k_index; k<nZ_d-1; k+=k_stride) {
      int i=calcIndex(n,k);
      dVardt[i] += -1*dfdz(var,0,k)*(n*M_PI/aspectRatio_d) * psi[i];
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
  for(int k=1+k_index; k<nZ_d-1; k+=k_stride) {
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
  for(int n=1+n_index; n<nN_d; n+=n_stride) {
    // Contribution FROM var[n>0] and vars.omg[n>0]
    int o;
    for(int m=1; m<n; ++m){
      // Case n = n' + n''
      o = n-m;
      for(int k=1+k_index; k<nZ_d-1; k+=k_stride) {
        int im = calcIndex(m,k);
        int io = calcIndex(o,k);
        int in = calcIndex(n,k);
        dVardt[in] +=
          -M_PI/(2.0*aspectRatio_d)*(
          -m*dfdz(psi,o,k)*var[im]
          +o*dfdz(var,m,k)*psi[io]
          );
      }
    }
    for(int m=n+1; m<nN_d; ++m){
      // Case n = n' - n''
      o = m-n;
      for(int k=1+k_index; k<nZ_d-1; k+=k_stride) {
        int im = calcIndex(m,k);
        int io = calcIndex(o,k);
        int in = calcIndex(n,k);
        dVardt[in] +=
          -M_PI/(2.0*aspectRatio_d)*(
          +m*dfdz(psi,o,k)*var[im]
          +o*dfdz(var,m,k)*psi[io]
          );
      }
    }
    for(int m=1; m+n<nN_d; ++m){
      // Case n= n'' - n'
      o = n+m;
      for(int k=1+k_index; k<nZ_d-1; k+=k_stride) {
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
            (var[calcIndex(i+1,k)]*-dfdz_s(psi,i+1,k) -
             var[calcIndex(i-1,k)]*-dfdz_s(psi,i-1,k))*oodx_d*0.5 +
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
{
  dt = c.initialDt;

  nonlinearTerm.initialiseData();
  nonlinearTerm.setupFFTW();

  thomasAlgorithm = new ThomasAlgorithmGPU(c.nZ, c.nN, c.aspectRatio, c.oodz2);

  gpuErrchk(cudaMemcpyToSymbol(nX_d, &c.nX, sizeof(c.nX), 0, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(nN_d, &c.nN, sizeof(c.nN), 0, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(nZ_d, &c.nZ, sizeof(c.nZ), 0, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(oodz_d, &c.oodz, sizeof(c.oodz), 0, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(oodx_d, &c.oodx, sizeof(c.oodx), 0, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(oodz2_d, &c.oodz2, sizeof(c.oodz2), 0, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(aspectRatio_d, &c.aspectRatio, sizeof(c.aspectRatio), 0, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(Ra_d, &c.Ra, sizeof(c.Ra), 0, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(Pr_d, &c.Pr, sizeof(c.Pr), 0, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(tau_d, &c.tau, sizeof(c.tau), 0, cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpyToSymbol(RaXi_d, &c.RaXi, sizeof(c.RaXi), 0, cudaMemcpyHostToDevice));
}

SimGPU::~SimGPU() {
  delete thomasAlgorithm;
}

void SimGPU::computeLinearTemperatureDerivative() {
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_computeLinearTemperatureDerivative<<<numBlocks,threadsPerBlock>>>(vars.dTmpdt.getCurrent(), vars.tmp.getCurrent());
}

void SimGPU::computeLinearVorticityDerivative() {
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_computeLinearVorticityDerivative<<<numBlocks,threadsPerBlock>>>(vars.dOmgdt.getCurrent(), vars.omg.getCurrent(), vars.tmp.getCurrent());
}

void SimGPU::computeLinearXiDerivative() {
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_computeLinearXiDerivative<<<numBlocks,threadsPerBlock>>>(vars.dXidt.getCurrent(), vars.xi.getCurrent(), vars.dOmgdt.getCurrent(), vars.omg.getCurrent());
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
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN - 1 + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  dim3 fillThreadsPerBlock(c.threadsPerBlock_x*c.threadsPerBlock_y);
  dim3 fillNumBlocks((c.nZ - 1 + fillThreadsPerBlock.x)/fillThreadsPerBlock.x);
  gpu_fillMode<<<fillNumBlocks,fillThreadsPerBlock>>>(vars.dOmgdt.getCurrent(), makeComplex(0.0, 0.0), 0);
  gpu_fillMode<<<fillNumBlocks,fillThreadsPerBlock>>>(vars.dTmpdt.getCurrent(), makeComplex(0.0, 0.0), 0);
  gpu_addAdvectionApproximation<<<numBlocks,threadsPerBlock>>>(
      vars.dTmpdt.getCurrent(), vars.tmp.getCurrent(), vars.psi.getCurrent());
  if(c.isDoubleDiffusion) {
    gpu_fillMode<<<fillNumBlocks,fillThreadsPerBlock>>>(vars.dXidt.getCurrent(), makeComplex(0.0, 0.0), 0);
    gpu_addAdvectionApproximation<<<numBlocks,threadsPerBlock>>>(
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
  vars.advanceDerivatives();
  solveForPsi();
}

void SimGPU::computeNonlinearDerivativeSpectralTransform(VariableGPU& dVardt, const VariableGPU& var) {
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nX - 1 + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_computeNonlinearDerivativeSpectralTransform<<<numBlocks,threadsPerBlock>>>(nonlinearTerm.spatialData_d, var.spatialData_d, vars.psi.spatialData_d);
  nonlinearTerm.toSpectral();
  gpu_plusEquals<<<numBlocks,threadsPerBlock>>>(dVardt.data_d, nonlinearTerm.data_d);
}

void SimGPU::computeNonlinearTemperatureDerivative() {
  // Calculate n=0 gpu_mode
  dim3 n0ThreadsPerBlock(c.threadsPerBlock_x*c.threadsPerBlock_y);
  dim3 n0NumBlocks((c.nZ - 1 + n0ThreadsPerBlock.x)/n0ThreadsPerBlock.x);
  gpu_computeNonlinearDerivativeN0<<<n0NumBlocks,n0ThreadsPerBlock>>>(vars.dTmpdt.getCurrent(), vars.tmp.getCurrent(), vars.psi.getCurrent());

  // Calculate other gpu_modes
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_computeNonlinearDerivative<<<numBlocks,threadsPerBlock>>>(vars.dTmpdt.getCurrent(), vars.tmp.getCurrent(), vars.psi.getCurrent(), -1);
}

void SimGPU::computeNonlinearXiDerivative() {
  // Calculate n=0 gpu_mode
  dim3 n0ThreadsPerBlock(c.threadsPerBlock_x*c.threadsPerBlock_y);
  dim3 n0NumBlocks((c.nZ - 1 + n0ThreadsPerBlock.x)/n0ThreadsPerBlock.x);
  gpu_computeNonlinearDerivativeN0<<<n0NumBlocks,n0ThreadsPerBlock>>>(vars.dXidt.getCurrent(), vars.xi.getCurrent(), vars.psi.getCurrent());

  // Calculate other gpu_modes
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_computeNonlinearDerivative<<<numBlocks,threadsPerBlock>>>(vars.dXidt.getCurrent(), vars.xi.getCurrent(), vars.psi.getCurrent(), -1);
}

void SimGPU::computeNonlinearVorticityDerivative() {
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN - 1 + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_computeNonlinearDerivative<<<numBlocks,threadsPerBlock>>>(vars.dOmgdt.getCurrent(), vars.omg.getCurrent(), vars.psi.getCurrent(), 1);
}

void SimGPU::runNonLinearStep(real f) {
  computeLinearDerivatives();
  computeNonlinearDerivatives();
  vars.updateVars(dt, f);
  vars.advanceDerivatives();
  solveForPsi();
}

//void SimGPU::computeNonlinearDerivatives() {
  //computeNonlinearTemperatureDerivative();
  //computeNonlinearVorticityDerivative();
  //if(c.isDoubleDiffusion) {
    //computeNonlinearXiDerivative();
  //}
//}

void SimGPU::computeNonlinearDerivatives() {
  vars.psi.toPhysical();

  gpuErrchk(cudaMemcpy(vars.psi.spatialData, vars.psi.spatialData_d, vars.psi.totalSize()*sizeof(vars.psi.spatialData[0]), cudaMemcpyDeviceToHost));
  for(int i=0; i<400; ++i) {
    std::cout << vars.psi.spatialData[i] << std::endl;
  }
  vars.tmp.toPhysical();
  gpuErrchk(cudaMemcpy(vars.tmp.spatialData, vars.tmp.spatialData_d, vars.tmp.totalSize()*sizeof(vars.tmp.spatialData[0]), cudaMemcpyDeviceToHost));
  for(int i=0; i<400; ++i) {
    std::cout << vars.tmp.spatialData[i] << std::endl;
  }
  vars.omg.toPhysical();
  gpuErrchk(cudaMemcpy(vars.omg.spatialData, vars.omg.spatialData_d, vars.omg.totalSize()*sizeof(vars.omg.spatialData[0]), cudaMemcpyDeviceToHost));
  for(int i=0; i<400; ++i) {
    std::cout << vars.omg.spatialData[i] << std::endl;
  }
  computeNonlinearDerivativeSpectralTransform(vars.dOmgdt, vars.omg);
  gpuErrchk(cudaMemcpy(vars.dOmgdt.data, vars.dOmgdt.data_d, vars.dOmgdt.totalSize()*sizeof(vars.dOmgdt.data[0]), cudaMemcpyDeviceToHost));
  for(int i=0; i<400; ++i) {
    std::cout << vars.dOmgdt.data[i] << std::endl;
  }
  computeNonlinearDerivativeSpectralTransform(vars.dTmpdt, vars.tmp);
  gpuErrchk(cudaMemcpy(vars.dTmpdt.data, vars.dTmpdt.data_d, vars.dTmpdt.totalSize()*sizeof(vars.dTmpdt.data[0]), cudaMemcpyDeviceToHost));
  for(int i=0; i<400; ++i) {
    std::cout << vars.dTmpdt.data[i] << std::endl;
  }
  //exit(0);
  //vars.dTmpdt.copyToHost();
  //for(int i=300; i<400; ++i) {
    //std::cout << vars.dTmpdt.data[i] << std::endl;
  //}
  if(c.isDoubleDiffusion) {
    vars.xi.toPhysical();
    computeNonlinearDerivativeSpectralTransform(vars.dXidt, vars.xi);
  }
}

void SimGPU::runNonLinear() {
  // Load initial conditions
  vars.load(c.icFile);

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
