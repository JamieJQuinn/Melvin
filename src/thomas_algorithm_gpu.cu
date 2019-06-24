#include "thomas_algorithm_gpu.hpp"

__global__
void solveThomasAlgorithm(real *sol, const real *rhs, const real *wk1, const real *wk2, const real *sub, const int n, const int nZ) {
  int iN = n*nZ;

  // Forward Subsitution
  sol[0] = rhs[0]*wk1[0+iN];
  for (int i=1; i<nZ; ++i) {
    sol[i] = (rhs[i] - sub[i-1]*sol[i-1])*wk1[i+iN];
  }
  // Backward Substitution
  for (int i=nZ-2; i>=0; --i) {
    sol[i] -= wk2[i+iN]*sol[i+1];
  }
}

void ThomasAlgorithmGPU::solve(real *sol, const real *rhs, const int n) const {
  solveThomasAlgorithm<<<1,1>>>((real*)sol, (real*)rhs, (real*)wk1, (real*)wk2, (real*)sub, n, nZ);
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

ThomasAlgorithmGPU::ThomasAlgorithmGPU(const int nZ, const int nN, const int a, const real oodz2):
  nZ(nZ),
  oodz2(oodz2)
  {
  cudaMallocManaged(&wk1, nZ*nN*sizeof(real));
  cudaMallocManaged(&wk2, nZ*nN*sizeof(real));
  cudaMallocManaged(&sub, nZ*sizeof(real));

  // Precalculate tridiagonal arrays
  real * dia = new real [nZ];
  real * sup = new real [nZ];
  for(int k=0; k<nZ; ++k) {
    sub[k] = sup[k] = -oodz2;
  }
  for(int n=0; n<nN; ++n) {
    for(int k=0; k<nZ; ++k){
      dia[k] = pow(M_PI/a*n, 2) + 2*oodz2;
    }
    dia[0] = dia[nZ-1] = 1.0;
    sub[nZ-2] = sup[0] = 0.0;
    formTriDiagonalArraysForN(
    sub, dia, sup,
    wk1+n*nZ, wk2+n*nZ);
  }

  delete [] dia;
  delete [] sup;
}

ThomasAlgorithmGPU::~ThomasAlgorithmGPU() {
  cudaFree(wk1);
  cudaFree(wk2);
  cudaFree(sub);
}
