#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include <precision.hpp>
#include <constants.hpp>
#include <variable.hpp>
#include <boundary_conditions.hpp>
#include <sim.hpp>

#include <iostream>
#include <cmath>

using namespace std;

void check_equal(const real x, const real y) {
  CHECK(x == Approx(y).margin(1e-13));
}

void require_equal(const real x, const real y) {
  REQUIRE(x == Approx(y).margin(1e-13));
}

void require_within_error(const real x, const real y) {
  REQUIRE(x == Approx(y).margin(1e-3));
}

TEST_CASE( "Test periodic boundary conditions", "[]" ) {
  Constants c;
  c.nN = 5;
  c.nZ = 10;
  c.aspectRatio = 1;
  c.verticalBoundaryConditions=BoundaryConditions::periodic;
  c.calculateDerivedConstants();

  Variable tmp(c);
  tmp.initialiseData();

  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      real x = k/(c.nZ-1);
      tmp(n,k) = (n+1);
    }
  }

  tmp.applyBoundaryConditions();

  for(int n=0; n<c.nN; ++n) {
    REQUIRE(tmp.dfdz(n, 0) == Approx(0.0));
    REQUIRE(tmp.dfdz(n, c.nZ-1) == Approx(0.0));
  }
}

TEST_CASE("Test discrete sine transform inverts", "[]") {
  Constants c("test_constants.json");

  Variable var(c);
  var.initialiseData();

  for(int ix=0; ix<var.nX; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix*c.aspectRatio/(var.nX-1.0);
      var.spatialData[var.calcIndex(ix, k)] = 0.0;
      for(int n=0; n<c.nN; ++n) {
        var.spatialData[var.calcIndex(ix, k)] += (k+1)*sin(M_PI * n * x/c.aspectRatio);
      }
    }
  }

  var.toSpectral(true);
  var.toPhysical(true);

  for(int ix=1; ix<var.nX-1; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix*c.aspectRatio/(var.nX-1.0);
      real sum = 0.0;
      for(int n=0; n<c.nN; ++n) {
        sum += (k+1)*sin(M_PI * n * x/c.aspectRatio);
      }
      //cout << ix << " " << k << endl;
      require_equal(var.spatialData[var.calcIndex(ix,k)], sum);
    }
  }
}

TEST_CASE("Test forward discrete sine transform", "[]") {
  Constants c("test_constants.json");

  Variable var(c);
  var.initialiseData();

  for(int ix=0; ix<var.nX; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix*c.aspectRatio/(var.nX-1.0);
      var.spatialData[var.calcIndex(ix, k)] = 0.0;
      for(int n=0; n<c.nN; ++n) {
        var.spatialData[var.calcIndex(ix, k)] += (k+1)*sin(M_PI * n * x/c.aspectRatio);
      }
    }
  }

  var.toSpectral(true);

  for(int k=0; k<c.nZ; ++k) {
    for(int n=1; n<c.nN; ++n) {
      require_equal(var(n, k), k+1);
    }
  }
}

TEST_CASE("Test backward discrete sine transform", "[]") {
  Constants c("test_constants.json");

  Variable var(c);
  var.initialiseData();

  for(int k=0; k<c.nZ; ++k) {
    real z = k/(c.nZ-1.0);
    for(int n=1; n<c.nN; ++n) {
      var(n,k) = n*sin(M_PI*z);
    }
  }

  var.toPhysical(true);

  for(int ix=1; ix<var.nX-1; ++ix) {
    real x = ix/(var.nX-1.0);
    for(int k=0; k<c.nZ; ++k) {
      real z = k/(c.nZ-1.0);
      real sum = 0.0;
      for(int n=0; n<c.nN; ++n) {
        sum += n*sin(M_PI*z)*sin(M_PI * n * x);
      }
      //cout << ix << " " << k << endl;
      //cout << var.spatialData[var.calcIndex(ix,k)] << " " << sum << endl;
      require_equal(var.spatialData[var.calcIndex(ix,k)], sum);
    }
  }

  for(int k=0; k<c.nZ; ++k) {
    require_equal(var.spatialData[var.calcIndex(0,k)], 0.0);
    require_equal(var.spatialData[var.calcIndex(var.nX-1,k)], 0.0);
  }
}

TEST_CASE("Test discrete cosine transform inverts", "[]") {
  Constants c("test_constants.json");

  Variable var(c);
  var.initialiseData();

  for(int ix=0; ix<var.nX; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix*c.aspectRatio/(var.nX-1.0);
      var.spatialData[var.calcIndex(ix, k)] = cos(M_PI * (k+1) * x/c.aspectRatio) + 2.0*cos(M_PI*(k+2)*x/c.aspectRatio);
    }
  }

  var.toSpectral(false);
  var.toPhysical(false);

  for(int ix=0; ix<var.nX; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix*c.aspectRatio/(var.nX-1.0);
      require_equal(var.spatialData[var.calcIndex(ix, k)], cos(M_PI * (k+1) * x/c.aspectRatio) + 2.0*cos(M_PI*(k+2)*x/c.aspectRatio));
    }
  }
}

TEST_CASE("Test forward discrete cosine transform", "[]") {
  Constants c("test_constants.json");

  Variable var(c);
  var.initialiseData();

  for(int ix=0; ix<var.nX; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix*c.aspectRatio/(var.nX-1.0);
      var.spatialData[var.calcIndex(ix, k)] = 5.0 + cos(M_PI * (k+1) * x/c.aspectRatio) + 2.0*cos(M_PI*(k+2)*x/c.aspectRatio);
    }
  }

  var.toSpectral(false);

  for(int k=0; k<c.nZ; ++k) {
    require_equal(var(0, k), 5.0);
    for(int n=1; n<k+1; ++n) {
      check_equal(var(n, k), 0.0);
    }

    require_equal(var((k+1), k), 1.0);
    require_equal(var((k+2), k), 2.0);

    for(int n=k+3; n<c.nN; ++n) {
      require_equal(var(n, k), 0.0);
    }
  }
}

TEST_CASE("Test backward discrete cosine transform", "[]") {
  Constants c("test_constants.json");

  Variable var(c);
  var.initialiseData();

  for(int k=0; k<c.nZ; ++k) {
    var(0, k) = 5.0;
    var(2, k) = 7.0;
  }

  var.toPhysical(false);

  for(int ix=0; ix<var.nX; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix/(var.nX-1.0);
      require_equal(var.spatialData[var.calcIndex(ix, k)], 5.0 + 7.0*cos(M_PI*2*x));
    }
  }
}

TEST_CASE("Test simple nonlinear multiplication and transform", "[]") {
  Constants c("test_constants.json");

  Variable var1(c);
  var1.initialiseData();

  for(int k=0; k<c.nZ; ++k) {
    var1(0, k) = 5.0;
    var1(2, k) = 2.0;
  }

  Variable var2(c);
  var2.initialiseData();

  for(int k=0; k<c.nZ; ++k) {
    var2(0, k) = 2.0;
  }

  var1.toPhysical(false);
  var2.toPhysical(false);

  Variable result(c);
  result.initialiseData();

  for(int k=0; k<c.nZ; ++k) {
    for(int i=0; i<result.nX; ++i) {
      result.spatialData[result.calcIndex(i,k)] = 
        var1.spatialData[var1.calcIndex(i,k)] *
        var2.spatialData[var2.calcIndex(i,k)];
    }
  }

  result.toSpectral(false);

  for(int k=0; k<c.nZ; ++k) {
    check_equal(result(0,k), 2.0*var1(0,k));
    check_equal(result(2,k), 2.0*var1(2,k));
  }
}

TEST_CASE("Test spectral transform calculation of sine-only nonlinear terms", "[]") {
  Constants c("test_constants.json");
  c.nN = 500;
  c.nZ = 1000;
  c.calculateDerivedConstants();

  Sim sim(c);
  sim.vars.load(c.icFile);

  int nX = sim.vars.dOmgdt.nX;

  for(int k=0; k<c.nZ; ++k) {
    real z = k/(c.nZ-1.0);
    for(int i=0; i<nX; ++i) {
      real x = i/(nX-1.0);
      sim.vars.psi.spatialData[sim.vars.psi.calcIndex(i,k)] = sin(M_PI*3.0*z)*sin(M_PI*x);
      sim.vars.omg.spatialData[sim.vars.omg.calcIndex(i,k)] = sin(M_PI*z)*sin(M_PI*x);
    }
  }

  sim.computeNonlinearVorticityDerivative();
  sim.vars.dOmgdt.toPhysical(true);

  Variable solution(c);
  solution.initialiseData();

  for(int i=1; i<nX-1; ++i) {
    real x = i/(nX-1.0);
    for(int k=1; k<c.nZ-1; ++k) {
      real z = k/(c.nZ-1.0);
          solution.spatialData[solution.calcIndex(i,k)] = 
            -0.5*M_PI*M_PI/c.aspectRatio*sin(2.0*M_PI*x)
            *(-3.0*cos(3.0*M_PI*z)*sin(M_PI*z)
              + sin(3.0*M_PI*z)*cos(M_PI*z));
    }
  }

  solution.toSpectral(true);

  for(int i=1; i<nX-1; ++i) {
    real x = i/(nX-1.0);
    for(int k=1; k<c.nZ-1; ++k) {
      real z = k/(c.nZ-1.0);
      //cout << i << " " << k << endl;
      real dpsidx = M_PI/c.aspectRatio*sin(3.0*M_PI*z)*cos(M_PI*x);
      real domgdx = M_PI/c.aspectRatio*sin(M_PI*z)*cos(M_PI*x);
      real dpsidz = 3.0*M_PI*cos(3.0*M_PI*z)*sin(M_PI*x);
      real domgdz = M_PI*cos(M_PI*z)*sin(M_PI*x);

      require_within_error(sim.vars.dOmgdt.spatialData[solution.calcIndex(i,k)],
          dpsidz*domgdx - dpsidx*domgdz);
    }
  }

  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      //cout << n << " " << k << endl;
      require_within_error(sim.vars.dOmgdt(n,k), solution(n,k));
    }
  }
}

TEST_CASE("Test spectral transform calculation of sine-only nonlinear term 1", "[]") {
  Constants c("test_constants.json");
  c.nN = 400;
  c.nZ = 1000;
  c.calculateDerivedConstants();

  Sim sim(c);
  sim.vars.load(c.icFile);

  int nX = sim.vars.dOmgdt.nX;

  for(int k=0; k<c.nZ; ++k) {
    real z = k/(c.nZ-1.0);
    for(int i=0; i<nX; ++i) {
      real x = i/(nX-1.0);
      sim.vars.psi.spatialData[sim.vars.psi.calcIndex(i,k)] = sin(M_PI*x);
      sim.vars.omg.spatialData[sim.vars.psi.calcIndex(i,k)] = sin(M_PI*z)*sin(M_PI*x);
    }
  }

  sim.computeNonlinearVorticityDerivative();
  sim.vars.dOmgdt.toPhysical(true);

  Variable solution(c);
  solution.initialiseData();

  for(int i=1; i<nX-1; ++i) {
    real x = i/(nX-1.0);
    for(int k=1; k<c.nZ-1; ++k) {
      real z = k/(c.nZ-1.0);
          solution.spatialData[solution.calcIndex(i,k)] = 
            -0.5*pow(M_PI,2)/c.aspectRatio*sin(2.0*M_PI*x)*cos(M_PI*z);
    }
  }

  solution.toSpectral(true);

  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      require_equal(sim.vars.dOmgdt(n,k), solution(n,k));
    }
  }
}

TEST_CASE("Test spectral transform calculation of sine-only nonlinear term 2", "[]") {
  Constants c("test_constants.json");
  c.nN = 1000;
  c.nZ = 2000;
  c.calculateDerivedConstants();

  Sim sim(c);
  sim.vars.load(c.icFile);

  int nX = sim.vars.dOmgdt.nX;

  for(int k=0; k<c.nZ; ++k) {
    real z = k/(c.nZ-1.0);
    for(int i=0; i<nX; ++i) {
      real x = i/(nX-1.0);
      sim.vars.omg.spatialData[sim.vars.omg.calcIndex(i,k)] = sin(M_PI*x);
      sim.vars.psi.spatialData[sim.vars.psi.calcIndex(i,k)] = sin(M_PI*z)*sin(M_PI*x);
    }
  }

  sim.computeNonlinearVorticityDerivative();
  sim.vars.dOmgdt.toPhysical(true);

  Variable solution(c);
  solution.initialiseData();

  for(int i=1; i<nX-1; ++i) {
    real x = i/(nX-1.0);
    for(int k=1; k<c.nZ-1; ++k) {
      real z = k/(c.nZ-1.0);
          solution.spatialData[solution.calcIndex(i,k)] = 
            0.5*pow(M_PI,2)/c.aspectRatio*sin(2.0*M_PI*x)*cos(M_PI*z);
    }
  }

  solution.toSpectral(true);

  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      require_equal(sim.vars.dOmgdt(n,k), solution(n,k));
    }
  }
}

TEST_CASE("Test spatial dfdz derivative", "[]") {
  Constants c("test_constants.json");
  c.nN = 200;
  c.nZ = 1000;
  c.calculateDerivedConstants();

  Variable var(c);
  var.initialiseData();

  for(int i=0; i<var.nX; ++i) {
    for(int k=0; k<c.nZ; ++k) {
      real z = k/(c.nZ-1.0);
      real x = i/(var.nX-1.0);
      var.spatialData[var.calcIndex(i,k)] = sin(3.0*M_PI*z)*sin(M_PI*x);
    }
  }

  for(int k=0; k<c.nZ; ++k) {
    real z = k/(c.nZ-1.0);
    var(1,k) = sin(3.0*M_PI*z);
  }

  for(int k=1; k<c.nZ-1; ++k) {
    real z = k/(c.nZ-1.0);
    require_within_error(var.dfdz(1,k), 3.0*M_PI*cos(3.0*M_PI*z));
  }

  for(int i=1; i<var.nX-1; ++i) {
    for(int k=1; k<c.nZ-1; ++k) {
      //cout << i << " " << k << endl;
      real z = k/(c.nZ-1.0);
      real x = i/(var.nX-1.0);
      require_within_error(var.dfdzSpatial(i,k), 3.0*M_PI*cos(3.0*M_PI*z)*sin(M_PI*x));
    }
  }
}

TEST_CASE("Test spatial dfdx derivative", "[]") {
  Constants c("test_constants.json");
  c.nN = 500;
  c.nZ = 50;
  c.calculateDerivedConstants();

  Variable var(c);
  var.initialiseData();

  for(int i=0; i<var.nX; ++i) {
    for(int k=0; k<c.nZ; ++k) {
      real z = k/(c.nZ-1.0);
      real x = i/(var.nX-1.0);
      var.spatialData[var.calcIndex(i,k)] = sin(M_PI*z)*sin(10.0*M_PI*x);
    }
  }

  for(int i=1; i<var.nX-1; ++i) {
    for(int k=1; k<c.nZ-1; ++k) {
      real z = k/(c.nZ-1.0);
      real x = i/(var.nX-1.0);
      require_within_error(var.dfdx(i,k), 10.0*M_PI/c.aspectRatio*sin(M_PI*z)*cos(10.0*M_PI*x));
    }
  }
}


void check_all_equal(Variable &var1, Variable &var2) {
  for(int n=0; n<var1.nN; ++n) {
    for(int k=1; k<var1.nZ-1; ++k) {
      cout << n << " " << k << endl;
      check_equal(var1(n,k), var2(n,k));
    }
  }
}

void require_all_equal(Variable &var1, Variable &var2) {
  for(int n=0; n<var1.nN; ++n) {
    for(int k=1; k<var1.nZ-1; ++k) {
      cout << n << " " << k << endl;
      require_equal(var1(n,k), var2(n,k));
    }
  }
}

TEST_CASE("Test reading and writing from file", "[]") {
  Constants c;
  c.nN = 5;
  c.nZ = 10;
  c.calculateDerivedConstants();

  Variable tmp(c);
  tmp.initialiseData();

  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      tmp(n,k) = n+1 + k;
    }
  }

  std::ofstream outfile ("./test.dat", std::ios::out | std::ios::binary);
  tmp.writeToFile(outfile);
  outfile.close();

  Variable tmp2(c);
  tmp2.initialiseData();

  std::ifstream infile ("./test.dat", std::ios::in | std::ios::binary);
  tmp2.readFromFile(infile);
  infile.close();
  
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      require_equal(tmp(n,k), tmp2(n,k));
    }
  }
}
