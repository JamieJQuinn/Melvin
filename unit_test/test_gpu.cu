#include "catch.hpp"
#include "thomas_algorithm_gpu.hpp"
#include "thomas_algorithm.hpp"

#include <iostream>

using namespace std;

TEST_CASE( "GPU loads and prints values correctly", "[gpu]" ) {
  int nZ = 5;

  real *rhsGPU;
  real *solGPU;

  cudaMallocManaged(&rhsGPU, nZ*sizeof(real));
  cudaMallocManaged(&solGPU, nZ*sizeof(real));

  for(int i=0; i<nZ; ++i) {
    rhsGPU[i] = i+1;
    solGPU[i] = 0.0;
  }

  cudaDeviceSynchronize();

  for(int i=0; i<nZ; ++i) {
    REQUIRE(rhsGPU[i] == real(i+1));
  }
}

TEST_CASE( "GPU Thomas algorithm solves a system correctly", "[gpu]" ) {
  int nZ = 5;
  int nN = 5;

  ThomasAlgorithm ta (nZ, nN, 1, 0.01f);
  ThomasAlgorithmGPU taGPU (nZ, nN, 1, 0.01f);

  real rhs [] = {1, 2, 3, 4, 5};
  real sol [nZ];

  ta.solve((real*)sol, (real*)rhs, 2);

  real *rhsGPU;
  real *solGPU;

  cudaMallocManaged(&rhsGPU, nZ*sizeof(real));
  cudaMallocManaged(&solGPU, nZ*sizeof(real));

  for(int i=0; i<nZ; ++i) {
    rhsGPU[i] = i+1;
    solGPU[i] = 0.0;
  }

  // Check precalculation works

  for(int i=0; i<nZ; ++i) {
    REQUIRE(ta.sub[i] == Approx(taGPU.sub[i]));
  }

  for(int n=0; n<nN; ++n) {
    for(int k=0; k<nZ; ++k) {
      REQUIRE(ta.wk1[k+n*nZ] == Approx(taGPU.wk1[k+n*nZ]));
      if(k < nZ-1){
        REQUIRE(ta.wk2[k+n*nZ] == Approx(taGPU.wk2[k+n*nZ]));
      }
    }
  }

  taGPU.solve((real*)solGPU, (real*)rhsGPU, 2);

  cudaDeviceSynchronize();

  for(int i=0; i<nZ; ++i) {
    REQUIRE(solGPU[i] == sol[i]);
  }
}
