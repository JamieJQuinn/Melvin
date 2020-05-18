#include <sim_gpu.hpp>

#include <precision.hpp>
#include <numerical_methods.hpp>
#include <complex_gpu.hpp>

#include <math.h>
#include <iostream>

using std::cout;
using std::endl;

__device__
gpu_mode dfdz2(const gpu_mode *data, const int n, const int k, const int nZ, const int oodz2) {
  return (data[calcIndex(n, k+1)] - 2.0f*data[calcIndex(n, k)] + data[calcIndex(n, k-1)])*oodz2;
}

__device__
gpu_mode dfdz(const gpu_mode *data, const int n, const int k, const int nZ, const int oodz) {
  return (data[calcIndex(n, k+1)] - data[calcIndex(n, k-1)])*oodz*0.5;
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
void gpu_computeLinearTemperatureDerivative(gpu_mode *dTmpdt, const gpu_mode *tmp,
    const int nN, const int nZ, const real aspectRatio, const real oodz2) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=n_index; n<nN; n+=n_stride) {
    for(int k=1+k_index; k<nZ-1; k+=k_stride) {
      int i=calcIndex(n,k);
      dTmpdt[i] = dfdz2(tmp, n, k, nZ, oodz2) - sqr(n*M_PI/aspectRatio)*tmp[i];
    }
  }
}

__global__
void gpu_computeLinearVorticityDerivative(gpu_mode *dOmgdt, const gpu_mode *omg, const gpu_mode *tmp,
    const int nN, const int nZ, const real aspectRatio, const real Ra, const real Pr, const real oodz2) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=n_index; n<nN; n+=n_stride) {
    for(int k=1+k_index; k<nZ-1; k+=k_stride) {
      int i=calcIndex(n,k);
      dOmgdt[i] =
        Pr*(dfdz2(omg,n,k,nZ,oodz2) - sqr(n*M_PI/aspectRatio)*omg[i])
        + Pr*Ra*(n*M_PI/aspectRatio)*tmp[i];
    }
  }
}

__global__
void gpu_computeLinearXiDerivative(gpu_mode *dXidt, const gpu_mode *xi, gpu_mode *dOmgdt, const gpu_mode *omg,
    const int nN, const int nZ, const real tau, const real aspectRatio, const real RaXi, const real Pr, const real oodz2) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=n_index; n<nN; n+=n_stride) {
    for(int k=1+k_index; k<nZ-1; k+=k_stride) {
      int i=calcIndex(n,k);
      dXidt[i] = tau*(dfdz2(xi, n, k, nZ, oodz2) - pow(n*M_PI/aspectRatio, 2)*xi[i]);
      dOmgdt[i] += -RaXi*tau*Pr*(n*M_PI/aspectRatio)*xi[i];
    }
  }
}

__global__
void gpu_fillMode(gpu_mode *data, const gpu_mode value, const int n, const int nZ) {
  int index = threadIdx.x;
  int stride = blockDim.x;
  for(int k=index; k<nZ; k+=stride) {
    int i=calcIndex(n,k);
    data[i] = value;
  }
}

__global__
void gpu_addAdvectionApproximation(
    gpu_mode *dVardt, const gpu_mode *var,
    const gpu_mode *psi,
    const int nN, const int nZ, const real aspectRatio, const real oodz) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=1+n_index; n<nN; n+=n_stride) {
    for(int k=1+k_index; k<nZ-1; k+=k_stride) {
      int i=calcIndex(n,k);
      dVardt[i] += -1*dfdz(var,0,k,nZ,oodz)*(n*M_PI/aspectRatio) * psi[i];
    }
  }
}

__global__
void gpu_computeNonlinearDerivativeN0(
    gpu_mode *dVardt, const gpu_mode *var,
    const gpu_mode *psi,
    const int nN, const int nZ, const real aspectRatio, const real oodz) {
  int k_index = blockIdx.x*blockDim.x + threadIdx.x;
  int k_stride = blockDim.x*gridDim.x;
  for(int k=1+k_index; k<nZ-1; k+=k_stride) {
    for(int n=1; n<nN; ++n) {
      // Contribution TO var[n=0]
      int i=calcIndex(n,k);
      dVardt[calcIndex(0,k)] +=
        -M_PI/(2*aspectRatio)*n*(
          dfdz(psi,n,k,nZ,oodz)*var[i] +
          dfdz(var,n,k,nZ,oodz)*psi[i]
          );
      // Contribution FROM var[n=0]
      dVardt[i] +=
        -n*M_PI/aspectRatio*psi[i]*dfdz(var,0,k,nZ,oodz);
    }
  }
}

__global__
void gpu_computeNonlinearDerivative(
    gpu_mode *dVardt, const gpu_mode *var,
    const gpu_mode *psi,
    const int nN, const int nZ, const real aspectRatio, const real oodz,
    const int vorticityFactor) {
  int n_index = blockIdx.x*blockDim.x + threadIdx.x;
  int n_stride = blockDim.x*gridDim.x;
  int k_index = blockIdx.y*blockDim.y + threadIdx.y;
  int k_stride = blockDim.y*gridDim.y;
  for(int n=1+n_index; n<nN; n+=n_stride) {
    // Contribution FROM var[n>0] and vars.omg[n>0]
    int o;
    for(int m=1; m<n; ++m){
      // Case n = n' + n''
      o = n-m;
      for(int k=1+k_index; k<nZ-1; k+=k_stride) {
        int im = calcIndex(m,k);
        int io = calcIndex(o,k);
        int in = calcIndex(n,k);
        dVardt[in] +=
          -M_PI/(2.0*aspectRatio)*(
          -m*dfdz(psi,o,k,nZ,oodz)*var[im]
          +o*dfdz(var,m,k,nZ,oodz)*psi[io]
          );
      }
    }
    for(int m=n+1; m<nN; ++m){
      // Case n = n' - n''
      o = m-n;
      for(int k=1+k_index; k<nZ-1; k+=k_stride) {
        int im = calcIndex(m,k);
        int io = calcIndex(o,k);
        int in = calcIndex(n,k);
        dVardt[in] +=
          -M_PI/(2.0*aspectRatio)*(
          +m*dfdz(psi,o,k,nZ,oodz)*var[im]
          +o*dfdz(var,m,k,nZ,oodz)*psi[io]
          );
      }
    }
    for(int m=1; m+n<nN; ++m){
      // Case n= n'' - n'
      o = n+m;
      for(int k=1+k_index; k<nZ-1; k+=k_stride) {
        int im = calcIndex(m,k);
        int io = calcIndex(o,k);
        int in = calcIndex(n,k);
        dVardt[in] +=
          vorticityFactor*M_PI/(2.0*aspectRatio)*(
          +m*dfdz(psi,o,k,nZ,oodz)*var[im]
          +o*dfdz(var,m,k,nZ,oodz)*psi[io]
          );
      }
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
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_computeLinearTemperatureDerivative<<<numBlocks,threadsPerBlock>>>(vars.dTmpdt.getCurrent(), vars.tmp.getCurrent(), c.nN, c.nZ, c.aspectRatio, c.oodz2);
}

void SimGPU::computeLinearVorticityDerivative() {
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_computeLinearVorticityDerivative<<<numBlocks,threadsPerBlock>>>(vars.dOmgdt.getCurrent(), vars.omg.getCurrent(), vars.tmp.getCurrent(),
    c.nN, c.nZ, c.aspectRatio, c.Ra, c.Pr, c.oodz2);
}

void SimGPU::computeLinearXiDerivative() {
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
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
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN - 1 + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  dim3 fillThreadsPerBlock(c.threadsPerBlock_x*c.threadsPerBlock_y);
  dim3 fillNumBlocks((c.nZ - 1 + fillThreadsPerBlock.x)/fillThreadsPerBlock.x);
  gpu_fillMode<<<fillNumBlocks,fillThreadsPerBlock>>>(vars.dOmgdt.getCurrent(), makeComplex(0.0, 0.0), 0, c.nZ);
  gpu_fillMode<<<fillNumBlocks,fillThreadsPerBlock>>>(vars.dTmpdt.getCurrent(), makeComplex(0.0, 0.0), 0, c.nZ);
  gpu_addAdvectionApproximation<<<numBlocks,threadsPerBlock>>>(
      vars.dTmpdt.getCurrent(), vars.tmp.getCurrent(), vars.psi.getCurrent(),
      c.nN, c.nZ, c.aspectRatio, c.oodz);
  if(c.isDoubleDiffusion) {
    gpu_fillMode<<<fillNumBlocks,fillThreadsPerBlock>>>(vars.dXidt.getCurrent(), makeComplex(0.0, 0.0), 0, c.nZ);
    gpu_addAdvectionApproximation<<<numBlocks,threadsPerBlock>>>(
        vars.dXidt.getCurrent(), vars.xi.getCurrent(), vars.psi.getCurrent(),
        c.nN, c.nZ, c.aspectRatio, c.oodz);
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

void SimGPU::computeNonlinearTemperatureDerivative() {
  // Calculate n=0 gpu_mode
  dim3 n0ThreadsPerBlock(c.threadsPerBlock_x*c.threadsPerBlock_y);
  dim3 n0NumBlocks((c.nZ - 1 + n0ThreadsPerBlock.x)/n0ThreadsPerBlock.x);
  gpu_computeNonlinearDerivativeN0<<<n0NumBlocks,n0ThreadsPerBlock>>>(vars.dTmpdt.getCurrent(), vars.tmp.getCurrent(), vars.psi.getCurrent(), 
      c.nN, c.nZ, c.aspectRatio, c.oodz);

  // Calculate other gpu_modes
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_computeNonlinearDerivative<<<numBlocks,threadsPerBlock>>>(vars.dTmpdt.getCurrent(), vars.tmp.getCurrent(), vars.psi.getCurrent(), 
      c.nN, c.nZ, c.aspectRatio, c.oodz, -1);
}

void SimGPU::computeNonlinearXiDerivative() {
  // Calculate n=0 gpu_mode
  dim3 n0ThreadsPerBlock(c.threadsPerBlock_x*c.threadsPerBlock_y);
  dim3 n0NumBlocks((c.nZ - 1 + n0ThreadsPerBlock.x)/n0ThreadsPerBlock.x);
  gpu_computeNonlinearDerivativeN0<<<n0NumBlocks,n0ThreadsPerBlock>>>(vars.dXidt.getCurrent(), vars.xi.getCurrent(), vars.psi.getCurrent(), 
      c.nN, c.nZ, c.aspectRatio, c.oodz);

  // Calculate other gpu_modes
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_computeNonlinearDerivative<<<numBlocks,threadsPerBlock>>>(vars.dXidt.getCurrent(), vars.xi.getCurrent(), vars.psi.getCurrent(), 
      c.nN, c.nZ, c.aspectRatio, c.oodz, -1);
}

void SimGPU::computeNonlinearVorticityDerivative() {
  dim3 threadsPerBlock(c.threadsPerBlock_x,c.threadsPerBlock_y);
  dim3 numBlocks((c.nN - 1 + threadsPerBlock.x - 1)/threadsPerBlock.x, (c.nZ - 2 + threadsPerBlock.y - 1)/threadsPerBlock.y);
  gpu_computeNonlinearDerivative<<<numBlocks,threadsPerBlock>>>(vars.dOmgdt.getCurrent(), vars.omg.getCurrent(), vars.psi.getCurrent(), 
      c.nN, c.nZ, c.aspectRatio, c.oodz, 1);
}

void SimGPU::runNonLinearStep(real f) {
  computeLinearDerivatives();
  computeNonlinearDerivatives();
  vars.updateVars(dt, f);
  vars.advanceDerivatives();
  solveForPsi();
}

void SimGPU::computeNonlinearDerivatives() {
  computeNonlinearTemperatureDerivative();
  computeNonlinearVorticityDerivative();
  if(c.isDoubleDiffusion) {
    computeNonlinearXiDerivative();
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
