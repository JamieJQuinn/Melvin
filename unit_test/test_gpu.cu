#include "catch.hpp"
#include "thomas_algorithm_gpu.hpp"
#include "thomas_algorithm.hpp"
#include "constants.hpp"
#include "variable.hpp"
#include "variable_gpu.hpp"
#include "variables.hpp"
#include "sim.hpp"
#include "sim_gpu.hpp"

#include <chrono>
#include <thread>
#include <iostream>

using namespace std;

using Clock = std::chrono::high_resolution_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using namespace std::literals::chrono_literals;

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
  int nZ = 10;
  int nN = 5;

  ThomasAlgorithm ta (nZ, nN, 1, 0.01f);
  ThomasAlgorithmGPU taGPU (nZ, nN, 1, 0.01f);

  real rhs [nZ];
  real sol [nZ];

  real *rhsGPU;
  real *solGPU;

  cudaMallocManaged(&rhsGPU, nN*nZ*sizeof(real));
  cudaMallocManaged(&solGPU, nN*nZ*sizeof(real));

  for(int i=0; i<nZ; ++i) {
    rhsGPU[i+2*nZ] = i+1;
    solGPU[i+2*nZ] = 0.0;
    rhs[i] = i+1;
    sol[i] = 0.0;
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

  ta.solve((real*)sol, (real*)rhs, 2);
  taGPU.solve((real*)solGPU, (real*)rhsGPU);

  cudaDeviceSynchronize();

  for(int i=0; i<nZ; ++i) {
    CHECK(solGPU[i+2*nZ] == Approx(sol[i]));
    CHECK(solGPU[i+0*nZ] == Approx(0.0));
    CHECK(solGPU[i+1*nZ] == Approx(0.0));
    CHECK(solGPU[i+3*nZ] == Approx(0.0));
    CHECK(solGPU[i+4*nZ] == Approx(0.0));
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

TEST_CASE("SimGPU initialises OK", "[gpu]") {
  Constants c;
  c.nN = 5;
  c.nZ = 10;
  c.aspectRatio = 1;
  c.calculateDerivedConstants();

  // Create GPU variables
  SimGPU sim(c);

  sim.vars.reinit(2.0);

  sim.vars.tmp(0,0) = 1.0;

  REQUIRE(sim.vars.tmp(0,0) == Approx(1.0));
  REQUIRE(sim.vars.tmp(0,1) == Approx(2.0));
}

TEST_CASE("Linear step calculates correctly", "[gpu]") {
  Constants c;
  c.nN = 5;
  c.nZ = 10;
  c.aspectRatio = 1.3;
  c.Pr = 1.0;
  c.Ra = 2.5;
  c.RaXi = 2.0;
  c.tau = 0.01;
  c.isDoubleDiffusion = true;
  c.calculateDerivedConstants();

  Sim s(c);
  SimGPU sGPU(c);

  // Load both with same test data
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      s.vars.omg(n,k) = (float)k;
      s.vars.tmp(n,k) = (float)k/c.nZ;
      s.vars.psi(n,k) = (float)k/c.nN;
      s.vars.xi(n,k) = (float)k/c.nN;
    }
  }

  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      sGPU.vars.omg(n,k) = (float)k;
      sGPU.vars.tmp(n,k) = (float)k/c.nZ;
      sGPU.vars.psi(n,k) = (float)k/c.nN;
      sGPU.vars.xi(n,k) = (float)k/c.nN;
    }
  }

  s.runLinearStep();
  sGPU.runLinearStep();
  cudaDeviceSynchronize();

  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      CHECK(sGPU.vars.dOmgdt(n,k) == Approx(s.vars.dOmgdt(n,k)));
      CHECK(sGPU.vars.dTmpdt(n,k) == Approx(s.vars.dTmpdt(n,k)));
      CHECK(sGPU.vars.dXidt(n,k) == Approx(s.vars.dXidt(n,k)));
      CHECK(sGPU.vars.tmp(n,k) == Approx(s.vars.tmp(n,k)));
      CHECK(sGPU.vars.psi(n,k) == Approx(s.vars.psi(n,k)));
      CHECK(sGPU.vars.omg(n,k) == Approx(s.vars.omg(n,k)));
      CHECK(sGPU.vars.xi(n,k) == Approx(s.vars.xi(n,k)));
    }
  }
}

TEST_CASE("Variable loads from file", "[gpu]") {
  Constants c;
  c.nN = 5;
  c.nZ = 10;
  c.aspectRatio = 1;
  c.calculateDerivedConstants();

  // Create GPU variables
  Variables<VariableGPU> varsGPU(c);
  Variables<Variable> vars(c);

  varsGPU.load("../test/benchmark/ICn1nZ101nN51");
  vars.load("../test/benchmark/ICn1nZ101nN51");

  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      CHECK(varsGPU.tmp(n, k) == Approx(vars.tmp(n, k)));
      CHECK(varsGPU.dTmpdt(n, k) == Approx(vars.dTmpdt(n, k)));
    }
  }
}

/*TEST_CASE("Rayleigh crit checker works with GPU", "[gpu]") {*/
  /*Constants c;*/
  /*c.nN = 5;*/
  /*c.nZ = 10;*/
  /*c.aspectRatio = 1.3;*/
  /*c.Pr = 1.0;*/
  /*c.Ra = 2.5;*/
  /*c.RaXi = 2.0;*/
  /*c.tau = 0.01;*/
  /*c.isDoubleDiffusion = true;*/
  /*c.isCudaEnabled = true;*/
  /*c.calculateDerivedConstants();*/

  /*CriticalRayleighChecker crc(c);*/

/*}*/

TEST_CASE("Benchmarking the linear step", "[gpu]") {
  Constants c;
  c.nN = 512;
  c.nZ = 1024;
  c.aspectRatio = 1.3;
  c.Pr = 1.0;
  c.Ra = 2.5;
  c.RaXi = 2.0;
  c.tau = 0.01;
  c.isDoubleDiffusion = true;
  c.calculateDerivedConstants();

  Sim s(c);

  c.isCudaEnabled = true;
  c.threadsPerBlock_x = 16;
  c.threadsPerBlock_y = 32;
  SimGPU sGPU(c);

  // Load both with same test data
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      s.vars.omg(n,k) = (float)k;
      s.vars.tmp(n,k) = (float)k/c.nZ;
      s.vars.psi(n,k) = (float)k/c.nN;
      s.vars.xi(n,k) = (float)k/c.nN;
    }
  }

  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      sGPU.vars.omg(n,k) = (float)k;
      sGPU.vars.tmp(n,k) = (float)k/c.nZ;
      sGPU.vars.psi(n,k) = (float)k/c.nN;
      sGPU.vars.xi(n,k) = (float)k/c.nN;
    }
  }

  time_point<Clock> start = Clock::now();
  s.computeLinearDerivatives();
  s.addAdvectionApproximation();
  s.vars.updateVars(s.dt);
  s.vars.advanceDerivatives();
  s.solveForPsi();
  time_point<Clock> end = Clock::now();
  std::chrono::duration<int64_t, std::nano> diff = end-start;
  cout << "CPU version of full linear step: " << diff.count() << endl;

  start = Clock::now();
  sGPU.computeLinearDerivatives();
  sGPU.addAdvectionApproximation();
  sGPU.vars.updateVars(sGPU.dt);
  sGPU.vars.advanceDerivatives();
  sGPU.solveForPsi();
  cudaDeviceSynchronize();
  end = Clock::now();
  diff = end-start;
  cout << "GPU version of full linear step: " << diff.count() << endl;

  start = Clock::now();
  s.computeLinearDerivatives();
  end = Clock::now();
  diff = end-start;
  cout << "CPU version of linear derivatives calculation: " << diff.count() << endl;

  start = Clock::now();
  sGPU.computeLinearDerivatives();
  cudaDeviceSynchronize();
  end = Clock::now();
  diff = end-start;
  cout << "GPU version of linear derivatives calculation: " << diff.count() << endl;

  start = Clock::now();
  s.addAdvectionApproximation();
  end = Clock::now();
  diff = end-start;
  cout << "CPU version of advection calculation: " << diff.count() << endl;

  start = Clock::now();
  sGPU.addAdvectionApproximation();
  cudaDeviceSynchronize();
  end = Clock::now();
  diff = end-start;
  cout << "GPU version of advection calculation: " << diff.count() << endl;

  start = Clock::now();
  s.vars.updateVars(s.dt);
  end = Clock::now();
  diff = end-start;
  cout << "CPU version of updating vars: " << diff.count() << endl;

  start = Clock::now();
  sGPU.vars.updateVars(sGPU.dt);
  cudaDeviceSynchronize();
  end = Clock::now();
  diff = end-start;
  cout << "GPU version of updating vars: " << diff.count() << endl;

  start = Clock::now();
  s.vars.advanceDerivatives();
  end = Clock::now();
  diff = end-start;
  cout << "CPU version of advancing derivatives: " << diff.count() << endl;

  start = Clock::now();
  sGPU.vars.advanceDerivatives();
  cudaDeviceSynchronize();
  end = Clock::now();
  diff = end-start;
  cout << "GPU version of advancing derivatives: " << diff.count() << endl;

  start = Clock::now();
  s.solveForPsi();
  end = Clock::now();
  diff = end-start;
  cout << "CPU version of Thomas algorithm: " << diff.count() << endl;

  start = Clock::now();
  sGPU.solveForPsi();
  cudaDeviceSynchronize();
  end = Clock::now();
  diff = end-start;
  cout << "GPU version of Thomas algorithm: " << diff.count() << endl;
}

