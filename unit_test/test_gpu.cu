#include "catch.hpp"
#include "thomas_algorithm_gpu.hpp"
#include "thomas_algorithm.hpp"
#include "constants.hpp"
#include "variable.hpp"
#include "variable_gpu.hpp"
#include "derivatives_gpu.hpp"
#include "derivatives.hpp"

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

TEST_CASE("GPU variable works", "[gpu]") {
  Constants c;
  c.nN = 5;
  c.nZ = 10;
  c.aspectRatio = 1;
  c.calculateDerivedConstants();

  // Create GPU variables
  VariableGPU tmp(c);
  tmp.initialiseData();

  tmp.data[5] = 1.0f;

  REQUIRE(tmp.data[5] == Approx(1.0f));
}

TEST_CASE("Temperature derivative calculates correctly", "[gpu]") {
  Constants c;
  c.nN = 5;
  c.nZ = 10;
  c.aspectRatio = 1;
  c.calculateDerivedConstants();

  // Create GPU variables
  VariableGPU tmpGPU(c);
  VariableGPU dTmpdtGPU(c, 2);
  tmpGPU.initialiseData();
  dTmpdtGPU.initialiseData();

  // Create CPU variables
  Variable tmp(c);
  Variable dTmpdt(c, 2);
  tmp.initialiseData();
  dTmpdt.initialiseData();

  // Load both with same test data
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      tmpGPU(n,k) = (float)k;
      tmp(n,k) = (float)k;
    }
  }

  computeLinearTemperatureDerivativeGPU(dTmpdtGPU, tmpGPU, c);
  computeLinearTemperatureDerivative(dTmpdt, tmp, c);

  cudaDeviceSynchronize();

  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      CHECK(dTmpdtGPU(n,k) == Approx(dTmpdt(n,k)));
    }
  }
}

TEST_CASE("Advection approximation works on the GPU", "[gpu]") {
  Constants c;
  c.nN = 5;
  c.nZ = 10;
  c.aspectRatio = 1;
  c.isDoubleDiffusion = true;
  c.calculateDerivedConstants();

  // Create GPU variables
  VariableGPU omgGPU(c);
  VariableGPU tmpGPU(c);
  VariableGPU psiGPU(c);
  VariableGPU dOmgdtGPU(c, 2);
  VariableGPU dTmpdtGPU(c, 2);
  VariableGPU xiGPU(c);
  VariableGPU dXidtGPU(c, 2);
  omgGPU.initialiseData();
  tmpGPU.initialiseData();
  psiGPU.initialiseData();
  dTmpdtGPU.initialiseData();
  dOmgdtGPU.initialiseData();
  xiGPU.initialiseData();
  dXidtGPU.initialiseData();

  // Create CPU variables
  Variable omg(c);
  Variable tmp(c);
  Variable psi(c);
  Variable dOmgdt(c, 2);
  Variable dTmpdt(c, 2);
  Variable xi(c);
  Variable dXidt(c, 2);
  omg.initialiseData();
  tmp.initialiseData();
  psi.initialiseData();
  dTmpdt.initialiseData();
  dOmgdt.initialiseData();
  xi.initialiseData();
  dXidt.initialiseData();

  // Load both with same test data
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      omgGPU(n,k) = (float)k;
      tmpGPU(n,k) = (float)k/c.nZ;
      psiGPU(n,k) = (float)k/c.nN;
      xiGPU(n,k) = (float)k/c.nN;
      omg(n,k) = (float)k;
      tmp(n,k) = (float)k/c.nZ;
      psi(n,k) = (float)k/c.nN;
      xi(n,k) = (float)k/c.nN;
    }
  }

  addAdvectionApproximationGPU(dTmpdtGPU, tmpGPU, dOmgdtGPU, omgGPU, dXidtGPU, xiGPU, psiGPU, c);
  addAdvectionApproximation(dTmpdt, tmp, dOmgdt, omg, dXidt, xi, psi, c);

  cudaDeviceSynchronize();

  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      CHECK(dTmpdtGPU(n,k) == Approx(dTmpdt(n,k)));
      CHECK(dXidtGPU(n,k) == Approx(dXidt(n,k)));
    }
  }
}

TEST_CASE("Vorticity derivative calculates correctly", "[gpu]") {
  Constants c;
  c.nN = 5;
  c.nZ = 10;
  c.aspectRatio = 1;
  c.calculateDerivedConstants();

  // Create GPU variables
  VariableGPU omgGPU(c);
  VariableGPU tmpGPU(c);
  VariableGPU dOmgdtGPU(c, 2);
  omgGPU.initialiseData();
  tmpGPU.initialiseData();
  dOmgdtGPU.initialiseData();

  // Create CPU variables
  Variable omg(c);
  Variable tmp(c);
  Variable dOmgdt(c, 2);
  omg.initialiseData();
  tmp.initialiseData();
  dOmgdt.initialiseData();

  // Load both with same test data
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      omgGPU(n,k) = (float)k;
      tmpGPU(n,k) = (float)k/c.nZ;
      omg(n,k) = (float)k;
      tmp(n,k) = (float)k/c.nZ;
    }
  }

  computeLinearVorticityDerivativeGPU(dOmgdtGPU, omgGPU, tmpGPU, c);
  computeLinearVorticityDerivative(dOmgdt, omg, tmp, c);

  cudaDeviceSynchronize();

  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      CHECK(dOmgdtGPU(n,k) == Approx(dOmgdt(n,k)));
    }
  }
}

TEST_CASE("Adams-Bashforth integrator works", "[gpu]") {
  Constants c;
  c.nN = 5;
  c.nZ = 10;
  c.aspectRatio = 1;
  c.calculateDerivedConstants();

  // Create GPU variables
  VariableGPU tmpGPU(c);
  VariableGPU dTmpdtGPU(c, 2);
  tmpGPU.initialiseData();
  dTmpdtGPU.initialiseData();

  // Create CPU variables
  Variable tmp(c);
  Variable dTmpdt(c, 2);
  tmp.initialiseData();
  dTmpdt.initialiseData();

  // Load both with same test data
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      tmpGPU(n,k) = (float)k;
      tmp(n,k) = (float)k;
    }
  }

  tmp.update(dTmpdt, 0.01, 1.0);
  tmpGPU.update(dTmpdtGPU, 0.01, 1.0);

  cudaDeviceSynchronize();

  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      CHECK(tmpGPU(n,k) == Approx(tmp(n,k)));
    }
  }
}
