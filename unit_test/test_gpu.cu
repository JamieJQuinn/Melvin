#include "catch.hpp"
#include "thomas_algorithm_gpu.hpp"
#include "thomas_algorithm.hpp"
#include "constants.hpp"
#include "variable.hpp"
#include "variable_gpu.hpp"
#include "variables.hpp"
#include "sim.hpp"
#include "sim_gpu.hpp"
#include "test_helper_functions.hpp"

#include <chrono>
#include <thread>
#include <iostream>

using Clock = std::chrono::high_resolution_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using namespace std::literals::chrono_literals;

using std::cout;
using std::endl;

using namespace std::complex_literals;


void test_main_vars(const Constants &c, const SimGPU &sGPU, const Sim &s, bool print=false) {
  for(auto variable : sGPU.vars.variableList) {
    variable->copyToHost();
  }
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      if(print) {
        cout << n << " " << k << endl;
      }
      require_equal(sGPU.vars.tmp(n,k), s.vars.tmp(n,k));
      require_equal(sGPU.vars.psi(n,k), s.vars.psi(n,k));
      require_equal(sGPU.vars.omg(n,k), s.vars.omg(n,k));
    }
  }
  if(c.isDoubleDiffusion) {
    for(int n=0; n<c.nN; ++n) {
      for(int k=0; k<c.nZ; ++k) {
        if(print) {
          cout << n << " " << k << endl;
        }
        require_equal(sGPU.vars.xi(n,k), s.vars.xi(n,k));
      }
    }
  }
}

void test_derivatives(const Constants &c, const SimGPU &sGPU, const Sim &s, bool print=false) {
  for(auto variable : sGPU.vars.variableList) {
    variable->copyToHost();
  }
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      if(print) {
        cout << n << " " << k << endl;
      }
      require_equal(sGPU.vars.dOmgdt(n,k), s.vars.dOmgdt(n,k));
      require_equal(sGPU.vars.dTmpdt(n,k), s.vars.dTmpdt(n,k));
    }
  }
  if(c.isDoubleDiffusion) {
    for(int n=0; n<c.nN; ++n) {
      for(int k=0; k<c.nZ; ++k) {
        if(print) {
          cout << n << " " << k << endl;
        }
        require_equal(sGPU.vars.dXidt(n,k), s.vars.dXidt(n,k));
      }
    }
  }
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      if(print) {
        cout << n << " " << k << endl;
      }
      require_equal(sGPU.vars.dOmgdt.getPrev(n,k), s.vars.dOmgdt.getPrev(n,k));
      require_equal(sGPU.vars.dTmpdt.getPrev(n,k), s.vars.dTmpdt.getPrev(n,k));
    }
  }
  if(c.isDoubleDiffusion) {
    for(int n=0; n<c.nN; ++n) {
      for(int k=0; k<c.nZ; ++k) {
        if(print) {
          cout << n << " " << k << endl;
        }
        require_equal(sGPU.vars.dXidt.getPrev(n,k), s.vars.dXidt.getPrev(n,k));
      }
    }
  }
}

void test_all_vars(const Constants &c, const SimGPU &sGPU, const Sim &s, bool print=false) {
  test_main_vars(c, sGPU, s, print);
  test_derivatives(c, sGPU, s, print);
}

//TEST_CASE( "GPU Thomas algorithm solves a system correctly", "[gpu]" ) {
  //int nZ = 10;
  //int nN = 5;

  //ThomasAlgorithm ta (nZ, nN, 1, 0.01f);
  //ThomasAlgorithmGPU taGPU (nZ, nN, 1, 0.01f);

  //real rhs [nZ];
  //real sol [nZ];

  //real *rhsGPU;
  //real *solGPU;

  //cudaMallocManaged(&rhsGPU, nN*nZ*sizeof(real));
  //cudaMallocManaged(&solGPU, nN*nZ*sizeof(real));

  //for(int i=0; i<nZ; ++i) {
    //rhsGPU[i+2*nZ] = i+1;
    //solGPU[i+2*nZ] = 0.0;
    //rhs[i] = i+1;
    //sol[i] = 0.0;
  //}

  //// Check precalculation works

  //for(int i=0; i<nZ; ++i) {
    //REQUIRE(ta.sub[i] == Approx(taGPU.sub[i]));
  //}

  //for(int n=0; n<nN; ++n) {
    //for(int k=0; k<nZ; ++k) {
      //REQUIRE(ta.wk1[k+n*nZ] == Approx(taGPU.wk1[k+n*nZ]));
      //if(k < nZ-1){
        //REQUIRE(ta.wk2[k+n*nZ] == Approx(taGPU.wk2[k+n*nZ]));
      //}
    //}
  //}

  //ta.solve((real*)sol, (real*)rhs, 2);
  //taGPU.solve((real*)solGPU, (real*)rhsGPU);

  //cudaDeviceSynchronize();

  //for(int i=0; i<nZ; ++i) {
    //CHECK(solGPU[i+2*nZ] == Approx(sol[i]));
    //CHECK(solGPU[i+0*nZ] == Approx(0.0));
    //CHECK(solGPU[i+1*nZ] == Approx(0.0));
    //CHECK(solGPU[i+3*nZ] == Approx(0.0));
    //CHECK(solGPU[i+4*nZ] == Approx(0.0));
  //}
//}

TEST_CASE("Make sure copy to and from device works", "[gpu]") {
  Constants c;
  c.nN = 5;
  c.nZ = 10;
  c.aspectRatio = 1;
  c.calculateDerivedConstants();

  // Create GPU variables
  VariableGPU tmp(c);

  for(int i=0; i<tmp.totalSize(); ++i) {
    tmp.data[i] = i;
  }

  tmp.copyToDevice();

  for(int i=0; i<tmp.totalSize(); ++i) {
    tmp.data[i] = 0.0f;
  }

  tmp.copyToHost();

  for(int i=0; i<tmp.totalSize(); ++i) {
    require_equal(tmp.data[i], i);
  }
}

TEST_CASE("GPU Variables class load from file", "[gpu]") {
  Constants c("test_constants_periodic_gpu.json");

  Sim s(c);
  SimGPU sGPU(c);

  s.vars.load(c.icFile);
  sGPU.vars.load(c.icFile);

  Variable& tmp = s.vars.tmp;
  VariableGPU& tmpGPU = sGPU.vars.tmp;

  tmpGPU.copyToHost();

  for(int i=0; i<tmp.totalSize(); ++i) {
    require_equal(tmp.data[i], tmpGPU.data[i]);
  }
}

TEST_CASE("Test cuFFT works at all" "[gpu]") {
  int nx = 256;
  int nn = nx/2 + 1;

  real* spatial = new real[nx];
  mode* spectral = new mode[nn];

  // Initialise spatial data

  for(int ix=0; ix<nx; ++ix) {
    //cout << ix << " " << k << endl;
    real x = real(ix)/nx;
    spatial[ix] = 
        5.0 +
        2.0*7.0*cos(2.0*M_PI*2*x) - 2.0*2.0*sin(2.0*M_PI*2*x) +
        2.0*2.0*cos(2.0*M_PI*3*x) + 2.0*10.0*sin(2.0*M_PI*3*x);
  }

  // Allocate space on GPU & copy test data

  real* spatial_d = nullptr;
  gpu_mode* spectral_d = nullptr;

  cudaMalloc(&spatial_d, sizeof(real)*nx);
  cudaMalloc(&spectral_d, sizeof(gpu_mode)*nn);

  cudaMemcpy(spatial_d, spatial, sizeof(real)*nx, cudaMemcpyHostToDevice);

  // Setup plan, run FFT & copy back to host

  cufftResult result;
  cufftHandle plan;

  result = cufftPlan1d(&plan, nx, CUFFT_D2Z, 1);
  REQUIRE(result == CUFFT_SUCCESS);

  result = cufftExecD2Z(plan, (cufftDoubleReal*)spatial_d, (cufftDoubleComplex*)spectral_d);
  REQUIRE(result == CUFFT_SUCCESS);

  cudaMemcpy(spectral, spectral_d, sizeof(mode)*nn, cudaMemcpyDeviceToHost);

  // Check

  REQUIRE(spectral[0].real() == Approx(5.0*nx));
  REQUIRE(spectral[2].real() == Approx(7.0*nx));
  REQUIRE(spectral[2].imag() == Approx(2.0*nx));
  REQUIRE(spectral[3].real() == Approx(2.0*nx));
  REQUIRE(spectral[3].imag() == Approx(-10.0*nx));

  // Cleanup
  cufftDestroy(plan);
  cudaFree(spatial_d);
  cudaFree(spectral_d);

  delete [] spatial;
  delete [] spectral;
}

TEST_CASE("Test cuFFT many plan works" "[gpu]") {
  // Define base number of modes
  int nn = 170;
  // Number of gridpoints to stabalise nonlinear term
  int nx = 3*nn+1;
  // Number of modes cufft needs to work with given nx (unrelated to actual number of modes we're using)
  int cufft_nn = nx/2+1;
  // Number of z-levels
  int nz = 256;
  // Number of ghost points (even so that we have alignment)
  int nG = 2;

  int spatialRowSize = nx + 2*nG;
  int spectralRowSize = cufft_nn + 2*nG;

  int spatialVarSize = spatialRowSize*(nz+2*nG);
  int spectralVarSize = spectralRowSize*(nz+2*nG);

  real* spatial = new real[spatialVarSize];
  mode* spectral = new mode[spectralVarSize];

  // Initialise spatial data

  for(int k=0; k<nz; ++k) {
    for(int ix=0; ix<nx; ++ix) {
      real x = real(ix)/nx;
      real z = 1.5*(k+1);
      int index = (k+nG)*spatialRowSize + ix + nG;
      spatial[index] = 
          z +
          2.0*z*cos(2.0*M_PI*2*x) - 2.0*2.0*sin(2.0*M_PI*2*x) +
          2.0*2.0*cos(2.0*M_PI*3*x) + 2.0*z*sin(2.0*M_PI*3*x);
    }
  }

  // Allocate space on GPU & copy test data

  real* spatial_d = nullptr;
  gpu_mode* spectral_d = nullptr;

  cudaMalloc(&spatial_d, sizeof(real)*spatialVarSize);
  cudaMalloc(&spectral_d, sizeof(gpu_mode)*spectralVarSize);

  cudaMemcpy(spatial_d, spatial, sizeof(real)*spatialVarSize, cudaMemcpyHostToDevice);

  // Setup plan, run FFT & copy back to host

  cufftResult result;
  cufftHandle plan;

  int rank = 1;
  int n[] = {nx};
  int inembed[] = {spatialRowSize};
  int istride = 1;
  int idist = spatialRowSize;
  int onembed[] = {spectralRowSize};
  int ostride = 1;
  int odist = spectralRowSize;
  int batch = nz;
  cufftType type = CUFFT_D2Z;

  result = cufftPlanMany(&plan,
      rank, n,
      inembed, istride, idist,
      onembed, ostride, odist,
      type, batch);
  REQUIRE(result == CUFFT_SUCCESS);

  real *spatialOffset_d = spatial_d + (0+nG)*spatialRowSize + 0 + nG;
  gpu_mode *spectralOffset_d = spectral_d + (0+nG)*spectralRowSize + 0 + nG;

  result = cufftExecD2Z(plan, (cufftDoubleReal*)spatialOffset_d, (cufftDoubleComplex*)spectralOffset_d);
  REQUIRE(result == CUFFT_SUCCESS);


  if (cudaDeviceSynchronize() != cudaSuccess){
    fprintf(stderr, "Cuda error: Failed to synchronize\n");
    return;
  }

  cudaMemcpy(spectral, spectral_d, sizeof(mode)*spectralVarSize, cudaMemcpyDeviceToHost);

  // Check

  for(int k=0; k<nz; ++k) {
    //cout << k << endl;
    real z = 1.5*(k+1);
    int index = (k+nG)*spectralRowSize + 0 + nG;
    check_equal(spectral[index].real(), (z*nx));

    index = (k+nG)*spectralRowSize + 2 + nG;
    check_equal(spectral[index].real(), (z*nx));
    check_equal(spectral[index].imag(), (2.0*nx));

    index = (k+nG)*spectralRowSize + 3 + nG;
    check_equal(spectral[index].real(), (2.0*nx));
    check_equal(spectral[index].imag(), (-z*nx));
  }

  // Cleanup
  cufftDestroy(plan);
  cudaFree(spatial_d);
  cudaFree(spectral_d);

  delete [] spatial;
  delete [] spectral;
}

TEST_CASE("Test variableGPU forward fft transform" "[gpu]") {
  Constants c("test_constants_periodic_gpu.json");

  VariableGPU var(c);

  // Initialise spatial data
  for(int k=0; k<c.nZ; ++k) {
    for(int ix=0; ix<c.nX; ++ix) {
      real x = real(ix)/c.nX;
      real z = 1.5*(k+1);
      var.spatial(ix,k) = 
          z +
          2.0*z*cos(2.0*M_PI*2*x) - 2.0*2.0*sin(2.0*M_PI*2*x) +
          2.0*2.0*cos(2.0*M_PI*3*x) + 2.0*z*sin(2.0*M_PI*3*x);
    }
  }

  var.copyToDevice(true);

  var.toSpectral();

  var.copyToHost();

  // Check

  for(int k=0; k<c.nZ; ++k) {
    real z = 1.5*(k+1);
    require_equal(var(0,k).real(), z);
    require_equal(var(2,k).real(), z);
    require_equal(var(2,k).imag(), 2.0);
    require_equal(var(3,k).real(), 2.0);
    require_equal(var(3,k).imag(), -z);
  }
}

TEST_CASE("Test variableGPU backward fft transform" "[gpu]") {
  Constants c("test_constants_periodic_gpu.json");

  VariableGPU var(c);

  for(int k=0; k<c.nZ; ++k) {
    real z = 2.0*(k+1);
    var(0, k) = z;
    var(2, k) = z + 2.0i;
    var(3, k) = 2.0 - z*1.0i;
  }

  var.copyToDevice();

  var.toPhysical();

  var.copyToHost(true);

  for(int k=0; k<c.nZ; ++k) {
    for(int ix=0; ix<var.nX; ++ix) {
      //cout << ix << " " << k << endl;
      real x = real(ix)/var.nX;
      real z = 2.0*(k+1);
      require_equal(var.spatial(ix,k),
          z +
          2.0*z*cos(2.0*M_PI*2*x) - 2.0*2.0*sin(2.0*M_PI*2*x) +
          2.0*2.0*cos(2.0*M_PI*3*x) + 2.0*z*sin(2.0*M_PI*3*x));
    }
  }
}

TEST_CASE("Test cuFFT discrete Fourier transform", "[gpu]") {
  Constants c("test_constants_periodic_gpu.json");

  VariableGPU var(c);

  for(int k=0; k<c.nZ; ++k) {
    var(0, k) = 5.0 + k;
    var(2, k) = 7.0 + 2.0i;
    var(3, k) = 2.0 - 10.0i;
  }

  var.copyToDevice();

  var.toPhysical();

  var.copyToHost(true);

  for(int k=0; k<c.nZ; ++k) {
    for(int ix=0; ix<var.nX; ++ix) {
      //cout << ix << " " << k << endl;
      real x = real(ix)/var.nX;
      require_equal(var.spatial(ix,k),
          5.0 + k +
          2.0*7.0*cos(2.0*M_PI*2*x) - 2.0*2.0*sin(2.0*M_PI*2*x) +
          2.0*2.0*cos(2.0*M_PI*3*x) + 2.0*10.0*sin(2.0*M_PI*3*x));
    }
  }

  var.toSpectral();

  var.copyToHost();

  for(int k=0; k<c.nZ; ++k) {
    require_equal(var(0, k), 5.0 + k);
    require_equal(var(2, k), 7.0 + 2.0i);
    require_equal(var(3, k), 2.0 - 10.0i);
  }
}

TEST_CASE("Test GPU spatial nonlinear derivative", "[]") {
  Constants c("test_constants_periodic_gpu.json");

  SimGPU sim(c);

  for(int i=0; i<c.nX; ++i) {
    for(int k=0; k<c.nZ; ++k) {
      real z = k*c.dz;
      real x = i*c.dx;
      sim.vars.psi.spatial(i,k) = (pow(x,2)+2.0*x-3.0)*(pow(z,2) + 5.0*z + 7.0);
      sim.vars.omg.spatial(i,k) = (pow(x,2)/3.0+x-2.0)*(2.0*pow(z,2) + 3.0*z + 10.0);
    }
  }

  sim.vars.psi.copyToDevice(true);
  sim.vars.omg.copyToDevice(true);

  sim.computeNonlinearDerivativeSpectralTransform(sim.vars.dOmgdt, sim.vars.omg);

  sim.nonlinearTerm.copyToHost(true);

  for(int i=1; i<c.nX-1; ++i) {
    for(int k=1; k<c.nZ-1; ++k) {
      real z = k*c.dz;
      real x = i*c.dx;
      real dpsidx = (2.0*x + 2.0)*(pow(z,2) + 5.0*z + 7.0);
      real dpsidz = (pow(x,2)+2.0*x-3.0)*(2.0*z + 5.0);
      real domgdx = (2.0/3.0*x + 1.0)*(2.0*pow(z,2) + 3.0*z + 10.0);
      real domgdz = (pow(x,2)/3.0+x-2.0)*(4.0*z + 3.0);

      CHECK(sim.nonlinearTerm.spatial(i,k)
        == Approx(-(-dpsidz*domgdx + dpsidx*domgdz)).margin(0.5));
    }
  }
}

TEST_CASE("Linear step calculates correctly", "[gpu]") {
  Constants c("test_constants_periodic.json");
  Constants cGPU("test_constants_periodic_gpu.json");

  Sim s(c);
  SimGPU sGPU(cGPU);

  // Load both with same test data
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      s.vars.omg(n,k) = (real)k + real(n)/c.nN*1.0i;
      s.vars.tmp(n,k) = (real)k + real(n)/c.nN*1.0i;
      s.vars.psi(n,k) = (real)k/c.nN;
      s.vars.xi(n,k) = (real)k/c.nN;
    }
  }

  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      sGPU.vars.omg(n,k) = (real)k + real(n)/c.nN*1.0i;
      sGPU.vars.tmp(n,k) = (real)k + real(n)/c.nN*1.0i;
      sGPU.vars.psi(n,k) = (real)k/c.nN;
      sGPU.vars.xi(n,k) = (real)k/c.nN;
    }
  }

  sGPU.vars.omg.copyToDevice();
  sGPU.vars.tmp.copyToDevice();
  sGPU.vars.psi.copyToDevice();

  for(int i=0; i<10; ++i) {
    s.runLinearStep();
    sGPU.runLinearStep();
    cudaDeviceSynchronize();
  }
  test_main_vars(c, sGPU, s);
}

//TEST_CASE("Each stage of linear step calculates correctly", "[gpu]") {
  //Constants c("test_constants_periodic.json");
  //Constants cGPU("test_constants_periodic_gpu.json");

  //Sim s(c);
  //SimGPU sGPU(cGPU);

  //// Load both with same test data
  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //s.vars.omg(n,k) = (real)k + real(n)/c.nN*1.0i;
      //s.vars.tmp(n,k) = (real)k + (real)n/c.nN*1.0i;
      //s.vars.psi(n,k) = (real)k/c.nN;
      //s.vars.xi(n,k) = (real)k/c.nN;
    //}
  //}

  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //sGPU.vars.omg(n,k) = (real)k + real(n)/c.nN*1.0i;
      //sGPU.vars.tmp(n,k) = (real)k + (real)n/c.nN*1.0i;
      //sGPU.vars.psi(n,k) = (real)k/c.nN;
      //sGPU.vars.xi(n,k) = (real)k/c.nN;
    //}
  //}

  //sGPU.vars.omg.copyToDevice();
  //sGPU.vars.tmp.copyToDevice();
  //sGPU.vars.psi.copyToDevice();

  //for(int i=0; i<10; ++i) {
    //s.computeLinearDerivatives();
    //sGPU.computeLinearDerivatives();
    //cudaDeviceSynchronize();
    //test_derivatives(c, sGPU, s);

    //cout << "Linear derivatives fine" << endl;

    //s.addAdvectionApproximation();
    //sGPU.addAdvectionApproximation();
    //cudaDeviceSynchronize();
    //test_derivatives(c, sGPU, s);

    //cout << "Advection approximation fine" << endl;

    //real dt = s.dt;
    //s.vars.updateVars(dt);
    //sGPU.vars.updateVars(dt);
    //cudaDeviceSynchronize();
    //test_main_vars(c, sGPU, s);

    //cout << "updating variables fine" << endl;

    //s.applyTemperatureBoundaryConditions();
    //sGPU.vars.tmp.applyVerticalBoundaryConditions();
    //cudaDeviceSynchronize();
    //test_main_vars(c, sGPU, s);

    //s.applyVorticityBoundaryConditions();
    //sGPU.vars.omg.applyVerticalBoundaryConditions();
    //cudaDeviceSynchronize();
    //test_main_vars(c, sGPU, s);

    //cout << "temp and vorticity bcs fine" << endl;

    //s.vars.advanceDerivatives();
    //sGPU.vars.advanceDerivatives();
    //cudaDeviceSynchronize();
    //test_derivatives(c, sGPU, s);

    //cout << "Advancing derivatives fine" << endl;

    //s.solveForPsi();
    //sGPU.solveForPsi();
    //cudaDeviceSynchronize();
    //test_main_vars(c, sGPU, s);

    //cout << "Solving for psi fine" << endl;

    //s.applyPsiBoundaryConditions();
    //sGPU.vars.psi.applyVerticalBoundaryConditions();
    //cudaDeviceSynchronize();
    //test_main_vars(c, sGPU, s);

    //cout << "psi bcs fine" << endl;
  //}
//}

TEST_CASE("Linear derivatives calculate correctly", "[gpu]") {
  Constants c("test_constants_periodic.json");
  Constants cGPU("test_constants_periodic_gpu.json");

  Sim s(c);
  SimGPU sGPU(cGPU);

  // Load both with same test data
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      s.vars.omg(n,k) = (real)k + real(n)/c.nN*1.0i;
      s.vars.tmp(n,k) = (real)k + (real)n/c.nN*1.0i;
      s.vars.psi(n,k) = (real)k/c.nN;
      s.vars.xi(n,k) = (real)k/c.nN;
    }
  }

  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      sGPU.vars.omg(n,k) = (real)k + real(n)/c.nN*1.0i;
      sGPU.vars.tmp(n,k) = (real)k + (real)n/c.nN*1.0i;
      sGPU.vars.psi(n,k) = (real)k/c.nN;
      sGPU.vars.xi(n,k) = (real)k/c.nN;
    }
  }

  sGPU.vars.omg.copyToDevice();
  sGPU.vars.tmp.copyToDevice();
  sGPU.vars.psi.copyToDevice();

  s.computeLinearDerivatives();
  sGPU.computeLinearDerivatives();
  cudaDeviceSynchronize();

  test_derivatives(c, sGPU, s);
}

TEST_CASE("Nonlinear step calculates correctly", "[gpu]") {
  Constants c("test_constants_periodic.json");
  Constants cGPU("test_constants_periodic_gpu.json");

  Sim s(c);
  SimGPU sGPU(cGPU);

  s.vars.load(c.icFile);
  sGPU.vars.load(c.icFile);

  for(int i=0; i<2; ++i) {
    //cout << i << endl;

    //cout << "computing linear derivatives" << endl;
    s.computeLinearDerivatives();
    sGPU.computeLinearDerivatives();
    test_derivatives(c, sGPU, s);

    //cout << "computing nonlinear derivatives" << endl;
    s.computeNonlinearDerivatives();
    sGPU.computeNonlinearDerivatives();
    test_derivatives(c, sGPU, s);

    //cout << "updating vars" << endl;
    real f = 1.0;
    real dt = s.dt;
    s.vars.updateVars(dt, f);
    sGPU.vars.updateVars(dt, f);
    test_all_vars(c, sGPU, s);

    //cout << "applying temp boundary conditions" << endl;
    s.applyTemperatureBoundaryConditions();
    sGPU.vars.tmp.applyVerticalBoundaryConditions();
    test_all_vars(c, sGPU, s);

    //cout << "applying vorticity boundary conditions" << endl;
    s.applyVorticityBoundaryConditions();
    sGPU.vars.omg.applyVerticalBoundaryConditions();
    test_all_vars(c, sGPU, s);

    //cout << "advancing derivatives" << endl;
    s.vars.advanceDerivatives();
    sGPU.vars.advanceDerivatives();
    test_derivatives(c, sGPU, s);

    //cout << "solving for psi" << endl;
    s.solveForPsi();
    sGPU.solveForPsi();
    test_all_vars(c, sGPU, s);

    //cout << "applying psi boundary conditions" << endl;
    s.applyPsiBoundaryConditions();
    sGPU.vars.psi.applyVerticalBoundaryConditions();
    test_all_vars(c, sGPU, s);
  }

  test_all_vars(c, sGPU, s);
}

//TEST_CASE("Linear step calculates correctly", "[gpu]") {
  //Constants c("test_constants.json");

  //Sim s(c);
  //SimGPU sGPU(c);

  //s.vars.load(c.icFile);
  //sGPU.vars.load(c.icFile);

  //for(int i=0; i<10; ++i) {
    //s.runLinearStep();
    //sGPU.runLinearStep();
    //cudaDeviceSynchronize();
  //}

  //test_all_vars(c, sGPU, s);
//}

//TEST_CASE("Linear derivatives calculate correctly", "[gpu]") {
  //Constants c("test_constants.json");

  //Sim s(c);
  //SimGPU sGPU(c);

  //s.vars.load(c.icFile);
  //sGPU.vars.load(c.icFile);

  //s.computeLinearDerivatives();
  //sGPU.computeLinearDerivatives();
  //cudaDeviceSynchronize();

  //test_all_vars(c, sGPU, s);
//}

//TEST_CASE("Nonlinear step calculates correctly", "[gpu]") {
  //Constants cGPU("test_constants_periodic_gpu.json");
  //Constants c("test_constants_periodic.json");

  //Sim s(c);
  //SimGPU sGPU(cGPU);

  //s.vars.load(c.icFile);
  //sGPU.vars.load(c.icFile);

  //for(int i=0; i<10; ++i) {
    //s.runNonLinearStep();
    //sGPU.runNonLinearStep();
  //}

  //test_all_vars(c, sGPU, s);
//}

TEST_CASE("GPU vertical boundary conditions work", "[gpu]") {
  Constants c("test_constants_periodic_gpu.json");

  SimGPU s(c);

  s.vars.load(c.icFile);

  s.vars.tmp(0,0) = 5.0;
  s.vars.tmp(1,0) = 5.0;
  s.vars.tmp(0,c.nZ-1) = 3.0;

  s.vars.tmp.copyToDevice();

  s.vars.tmp.applyVerticalBoundaryConditions();

  s.vars.tmp.copyToHost();

  require_equal(s.vars.tmp(0,0), 1.0);
  require_equal(s.vars.tmp(1,0), 0.0);
  require_equal(s.vars.tmp(0,c.nZ-1), 0.0);
}

TEST_CASE("CPU and GPU spatial boundary conditions match", "[gpu]") {
  Constants c("test_constants_periodic.json");
  Constants cGPU("test_constants_periodic_gpu.json");

  Sim s(c);
  SimGPU sGPU(cGPU);

  // Load both with same test data
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      s.vars.omg(n,k) = (real)k + real(n)/c.nN*1.0i;
      s.vars.tmp(n,k) = (real)k + real(n)/c.nN*1.0i;
      s.vars.psi(n,k) = (real)k/c.nN;
      s.vars.xi(n,k) = (real)k/c.nN;
    }
  }

  s.vars.tmp.toPhysical();
  s.vars.omg.toPhysical();
  s.vars.psi.toPhysical();
  s.applyPhysicalBoundaryConditions();

  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      sGPU.vars.omg(n,k) = (real)k + real(n)/c.nN*1.0i;
      sGPU.vars.tmp(n,k) = (real)k + (real)n/c.nN*1.0i;
      sGPU.vars.psi(n,k) = (real)k/c.nN;
      sGPU.vars.xi(n,k) = (real)k/c.nN;
    }
  }

  sGPU.vars.omg.copyToDevice();
  sGPU.vars.tmp.copyToDevice();
  sGPU.vars.psi.copyToDevice();

  sGPU.vars.omg.toPhysical();
  sGPU.vars.tmp.toPhysical();
  sGPU.vars.psi.toPhysical();

  sGPU.vars.psi.applyPhysicalHorizontalBoundaryConditions();
  sGPU.vars.tmp.applyPhysicalHorizontalBoundaryConditions();
  sGPU.vars.omg.applyPhysicalHorizontalBoundaryConditions();

  sGPU.vars.omg.copyToHost(true);
  sGPU.vars.tmp.copyToHost(true);
  sGPU.vars.psi.copyToHost(true);

  for(int k=0; k<c.nZ; ++k) {
    for(int i=-1; i<c.nX+1; ++i) {
      require_equal(sGPU.vars.tmp.spatial(i,k), s.vars.tmp.spatial(i,k));
      require_equal(sGPU.vars.omg.spatial(i,k), s.vars.omg.spatial(i,k));
      require_equal(sGPU.vars.psi.spatial(i,k), s.vars.psi.spatial(i,k));
    }
  }
}

TEST_CASE("Nonlinear calculation matches pre-FFT", "[gpu]") {
  Constants c("test_constants_periodic.json");
  Constants cGPU("test_constants_periodic_gpu.json");

  Sim s(c);
  SimGPU sGPU(cGPU);

  // Load both with same test data
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      s.vars.omg(n,k) = (real)k + real(n)/c.nN*1.0i;
      s.vars.tmp(n,k) = (real)k + real(n)/c.nN*1.0i;
      s.vars.psi(n,k) = (real)k/c.nN;
      s.vars.xi(n,k) = (real)k/c.nN;
    }
  }

  s.vars.tmp.toPhysical();
  s.vars.omg.toPhysical();
  s.vars.psi.toPhysical();
  s.applyPhysicalBoundaryConditions();
  s.computeNonlinearVorticityDerivative();

  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      sGPU.vars.omg(n,k) = (real)k + real(n)/c.nN*1.0i;
      sGPU.vars.tmp(n,k) = (real)k + real(n)/c.nN*1.0i;
      sGPU.vars.psi(n,k) = (real)k/c.nN;
      sGPU.vars.xi(n,k) = (real)k/c.nN;
    }
  }

  sGPU.vars.omg.copyToDevice();
  sGPU.vars.tmp.copyToDevice();
  sGPU.vars.psi.copyToDevice();

  sGPU.vars.omg.toPhysical();
  sGPU.vars.tmp.toPhysical();
  sGPU.vars.psi.toPhysical();

  sGPU.vars.psi.applyPhysicalHorizontalBoundaryConditions();
  sGPU.vars.tmp.applyPhysicalHorizontalBoundaryConditions();
  sGPU.vars.omg.applyPhysicalHorizontalBoundaryConditions();
  sGPU.computeNonlinearDerivativeSpectralTransform(sGPU.vars.dOmgdt, sGPU.vars.omg);

  sGPU.nonlinearTerm.copyToHost(true);

  for(int k=0; k<c.nZ; ++k) {
    for(int i=0; i<c.nX; ++i) {
      //cout << k << " " << i << endl;
      require_equal(sGPU.nonlinearTerm.spatial(i,k), s.nonlinearSineTerm.spatial(i,k));
    }
  }
}

//TEST_CASE("Nonlinear temperature derivative calculates correctly", "[gpu]") {
  //Constants c("test_constants_periodic.json");
  //Constants cGPU("test_constants_periodic_gpu.json");

  //c.nN = 400;
  //cGPU.nN = c.nN;
  //c.nZ = 2*c.nN;
  //cGPU.nZ = c.nZ;
  //c.nX = 3*c.nN + 1;
  //cGPU.nX = c.nX;

  //c.calculateDerivedConstants();
  //cGPU.calculateDerivedConstants();

  //Sim s(c);
  //SimGPU sGPU(cGPU);

  //// Load both with same test data
  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //s.vars.omg(n,k) = (real)k;
      //s.vars.tmp(n,k) = (real)k + (real)n/c.nN*1.0i;
      //s.vars.psi(n,k) = (real)k/c.nN;
      //s.vars.xi(n,k) = (real)k/c.nN;
    //}
  //}

  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //sGPU.vars.omg(n,k) = (real)k;
      //sGPU.vars.tmp(n,k) = (real)k + (real)n/c.nN*1.0i;
      //sGPU.vars.psi(n,k) = (real)k/c.nN;
      //sGPU.vars.xi(n,k) = (real)k/c.nN;
    //}
  //}

  //sGPU.vars.omg.copyToDevice();
  //sGPU.vars.tmp.copyToDevice();
  //sGPU.vars.psi.copyToDevice();

  //time_point<Clock> start = Clock::now();
  //s.computeNonlinearDerivatives();
  //time_point<Clock> end = Clock::now();
  //std::chrono::duration<int64_t, std::nano> diff = end-start;
  //cout << "CPU version of nonlinear derivatives calculation: " << diff.count() << endl;

  //start = Clock::now();
  //sGPU.computeNonlinearDerivatives();
  //cudaDeviceSynchronize();
  //end = Clock::now();
  //diff = end-start;
  //cout << "GPU version of nonlinear derivatives calculation: " << diff.count() << endl;

  //sGPU.vars.dTmpdt.copyToHost();
  //sGPU.vars.dOmgdt.copyToHost();

  //for(int k=0; k<c.nZ; ++k) {
    //for(int n=0; n<c.nN; ++n) {
      //cout << k << " " << n << endl;
      //cout << sGPU.vars.dTmpdt(n,k) << s.vars.dTmpdt(n,k) << endl;
      //require_within_error(sGPU.vars.dTmpdt(n,k), s.vars.dTmpdt(n,k), 1e-10);
      //cout << sGPU.vars.dOmgdt(n,k) << s.vars.dOmgdt(n,k) << endl;
      //require_within_error(sGPU.vars.dOmgdt(n,k), s.vars.dOmgdt(n,k), 1e-10);
    //}
  //}
//}

//TEST_CASE("Nonlinear vorticity derivative calculates correctly", "[gpu]") {
  //Constants c;
  //c.nN = 64;
  //c.nZ = 128;
  //c.aspectRatio = 1.3;
  //c.Pr = 1.0;
  //c.Ra = 2.5;
  //c.RaXi = 2.0;
  //c.tau = 0.01;
  //c.isDoubleDiffusion = true;
  //c.calculateDerivedConstants();

  //Sim s(c);

  //c.isCudaEnabled = true;
  //c.threadsPerBlock_x = 16;
  //c.threadsPerBlock_y = 32;
  //SimGPU sGPU(c);

  //// Load both with same test data
  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //s.vars.omg(n,k) = (float)k;
      //s.vars.tmp(n,k) = (float)k;
      //s.vars.psi(n,k) = (float)k/c.nN;
      //s.vars.xi(n,k) = (float)k/c.nN;
    //}
  //}

  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //sGPU.vars.omg(n,k) = (float)k;
      //sGPU.vars.tmp(n,k) = (float)k;
      //sGPU.vars.psi(n,k) = (float)k/c.nN;
      //sGPU.vars.xi(n,k) = (float)k/c.nN;
    //}
  //}

  //time_point<Clock> start = Clock::now();
  //s.computeLinearVorticityDerivative();
  //s.computeNonlinearVorticityDerivative();
  //time_point<Clock> end = Clock::now();
  //std::chrono::duration<int64_t, std::nano> diff = end-start;
  //cout << "CPU version of nonlinear vorticity derivatives calculation: " << diff.count() << endl;

  //start = Clock::now();
  //sGPU.computeLinearVorticityDerivative();
  //sGPU.computeNonlinearVorticityDerivative();
  //cudaDeviceSynchronize();
  //end = Clock::now();
  //diff = end-start;
  //cout << "GPU version of nonlinear vorticity derivatives calculation: " << diff.count() << endl;

  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //REQUIRE(sGPU.vars.dOmgdt(n,k) == Approx(s.vars.dOmgdt(n,k)));
    //}
  //}
//}

//TEST_CASE("Nonlinear xi derivative calculates correctly", "[gpu]") {
  //Constants c;
  //c.nN = 64;
  //c.nZ = 128;
  //c.aspectRatio = 1.3;
  //c.Pr = 1.0;
  //c.Ra = 2.5;
  //c.RaXi = 2.0;
  //c.tau = 0.01;
  //c.isDoubleDiffusion = true;
  //c.calculateDerivedConstants();

  //Sim s(c);

  //c.isCudaEnabled = true;
  //c.threadsPerBlock_x = 16;
  //c.threadsPerBlock_y = 32;
  //SimGPU sGPU(c);

  //// Load both with same test data
  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //s.vars.omg(n,k) = (float)k;
      //s.vars.tmp(n,k) = (float)k/c.nZ;
      //s.vars.psi(n,k) = (float)k/c.nN;
      //s.vars.xi(n,k) = (float)k/c.nN;
    //}
  //}

  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //sGPU.vars.omg(n,k) = (float)k;
      //sGPU.vars.tmp(n,k) = (float)k/c.nZ;
      //sGPU.vars.psi(n,k) = (float)k/c.nN;
      //sGPU.vars.xi(n,k) = (float)k/c.nN;
    //}
  //}

  //time_point<Clock> start = Clock::now();
  //s.computeNonlinearXiDerivative();
  //time_point<Clock> end = Clock::now();
  //std::chrono::duration<int64_t, std::nano> diff = end-start;
  //cout << "CPU version of nonlinear xi derivatives calculation: " << diff.count() << endl;

  //start = Clock::now();
  //sGPU.computeNonlinearXiDerivative();
  //cudaDeviceSynchronize();
  //end = Clock::now();
  //diff = end-start;
  //cout << "GPU version of nonlinear xi derivatives calculation: " << diff.count() << endl;

  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //REQUIRE(sGPU.vars.dXidt(n,k) == Approx(s.vars.dXidt(n,k)));
    //}
  //}
//}

//TEST_CASE("Variable loads from file", "[gpu]") {
  //Constants c;
  //c.nN = 5;
  //c.nZ = 10;
  //c.aspectRatio = 1;
  //c.calculateDerivedConstants();

  //// Create GPU variables
  //Variables<VariableGPU> varsGPU(c);
  //Variables<Variable> vars(c);

  //varsGPU.load("../test/benchmark/ICn1nZ101nN51");
  //vars.load("../test/benchmark/ICn1nZ101nN51");

  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //CHECK(varsGPU.tmp(n, k) == Approx(vars.tmp(n, k)));
      //CHECK(varsGPU.dTmpdt(n, k) == Approx(vars.dTmpdt(n, k)));
    //}
  //}
//}

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

//TEST_CASE("Benchmarking the linear step", "[gpu]") {
  //Constants c;
  //c.nN = 512;
  //c.nZ = 1024;
  //c.aspectRatio = 1.3;
  //c.Pr = 1.0;
  //c.Ra = 2.5;
  //c.RaXi = 2.0;
  //c.tau = 0.01;
  //c.isDoubleDiffusion = true;
  //c.calculateDerivedConstants();

  //Sim s(c);

  //c.isCudaEnabled = true;
  //c.threadsPerBlock_x = 16;
  //c.threadsPerBlock_y = 32;
  //SimGPU sGPU(c);

  //// Load both with same test data
  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //s.vars.omg(n,k) = (float)k;
      //s.vars.tmp(n,k) = (float)k/c.nZ;
      //s.vars.psi(n,k) = (float)k/c.nN;
      //s.vars.xi(n,k) = (float)k/c.nN;
    //}
  //}

  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //sGPU.vars.omg(n,k) = (float)k;
      //sGPU.vars.tmp(n,k) = (float)k/c.nZ;
      //sGPU.vars.psi(n,k) = (float)k/c.nN;
      //sGPU.vars.xi(n,k) = (float)k/c.nN;
    //}
  //}

  //time_point<Clock> start = Clock::now();
  //s.computeLinearDerivatives();
  //s.addAdvectionApproximation();
  //s.vars.updateVars(s.dt);
  //s.vars.advanceDerivatives();
  //s.solveForPsi();
  //time_point<Clock> end = Clock::now();
  //std::chrono::duration<int64_t, std::nano> diff = end-start;
  //cout << "CPU version of full linear step: " << diff.count() << endl;

  //start = Clock::now();
  //sGPU.computeLinearDerivatives();
  //sGPU.addAdvectionApproximation();
  //sGPU.vars.updateVars(sGPU.dt);
  //sGPU.vars.advanceDerivatives();
  //sGPU.solveForPsi();
  //cudaDeviceSynchronize();
  //end = Clock::now();
  //diff = end-start;
  //cout << "GPU version of full linear step: " << diff.count() << endl;

  //start = Clock::now();
  //s.computeLinearDerivatives();
  //end = Clock::now();
  //diff = end-start;
  //cout << "CPU version of linear derivatives calculation: " << diff.count() << endl;

  //start = Clock::now();
  //sGPU.computeLinearDerivatives();
  //cudaDeviceSynchronize();
  //end = Clock::now();
  //diff = end-start;
  //cout << "GPU version of linear derivatives calculation: " << diff.count() << endl;

  //start = Clock::now();
  //s.addAdvectionApproximation();
  //end = Clock::now();
  //diff = end-start;
  //cout << "CPU version of advection calculation: " << diff.count() << endl;

  //start = Clock::now();
  //sGPU.addAdvectionApproximation();
  //cudaDeviceSynchronize();
  //end = Clock::now();
  //diff = end-start;
  //cout << "GPU version of advection calculation: " << diff.count() << endl;

  //start = Clock::now();
  //s.vars.updateVars(s.dt);
  //end = Clock::now();
  //diff = end-start;
  //cout << "CPU version of updating vars: " << diff.count() << endl;

  //start = Clock::now();
  //sGPU.vars.updateVars(sGPU.dt);
  //cudaDeviceSynchronize();
  //end = Clock::now();
  //diff = end-start;
  //cout << "GPU version of updating vars: " << diff.count() << endl;

  //start = Clock::now();
  //s.vars.advanceDerivatives();
  //end = Clock::now();
  //diff = end-start;
  //cout << "CPU version of advancing derivatives: " << diff.count() << endl;

  //start = Clock::now();
  //sGPU.vars.advanceDerivatives();
  //cudaDeviceSynchronize();
  //end = Clock::now();
  //diff = end-start;
  //cout << "GPU version of advancing derivatives: " << diff.count() << endl;

  //start = Clock::now();
  //s.solveForPsi();
  //end = Clock::now();
  //diff = end-start;
  //cout << "CPU version of Thomas algorithm: " << diff.count() << endl;

  //start = Clock::now();
  //sGPU.solveForPsi();
  //cudaDeviceSynchronize();
  //end = Clock::now();
  //diff = end-start;
  //cout << "GPU version of Thomas algorithm: " << diff.count() << endl;
//}

TEST_CASE("Test GPU complex poisson solver", "[]") {
  Constants c("test_constants_periodic_gpu.json");

  SimGPU sim(c);
  VariableGPU psi(c);

  for(int k=0; k<c.nZ; ++k) {
    for(int n=0; n<c.nN; ++n) {
      real z = c.dz*k;
      psi(n,k) = 2.0*sin(M_PI*z)*(n+1) + 3.0i*sin(2.0*M_PI*z);
    }
  }

  for(int n=0; n<c.nN; ++n) {
    sim.vars.omg(n,0) = sim.vars.omg(n,c.nZ-1) = 0.0;
    for(int k=1; k<c.nZ-1; ++k) {
      sim.vars.omg(n,k) = -(psi.dfdz2(n,k) - pow(real(n)*c.wavelength, 2)*psi(n,k));
    }
  }

  psi.copyToDevice();
  sim.vars.omg.copyToDevice();

  sim.solveForPsi();

  sim.vars.psi.copyToHost();

  for(int k=0; k<c.nZ; ++k) {
    for(int n=0; n<c.nN; ++n) {
      //cout << n << " " << k << endl;
      require_within_error(sim.vars.psi(n,k), psi(n,k));
    }
  }
}
