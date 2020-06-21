#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include <precision.hpp>
#include <constants.hpp>
#include <variable.hpp>
#include <boundary_conditions.hpp>
#include <sim.hpp>

#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using namespace std::complex_literals;

void require_within_error(const real x, const real y, const real margin=1e-3) {
  //cout << x << " " << y << endl;
  REQUIRE(x == Approx(y).margin(margin));
}

void require_within_error(const mode x, const mode y, const real margin=1e-3) {
  require_within_error(x.real(), y.real(), margin);
  require_within_error(x.imag(), y.imag(), margin);
}

void check_equal(const real x, const real y) {
  CHECK(x == Approx(y).margin(1e-13));
}

void require_equal(const real x, const real y) {
  REQUIRE(x == Approx(y).margin(1e-13));
}

void check_equal(const mode x, const mode y) {
  check_equal(x.real(), y.real());
  check_equal(x.imag(), y.imag());
}

void require_equal(const mode x, const mode y) {
  require_equal(x.real(), y.real());
  require_equal(x.imag(), y.imag());
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

TEST_CASE("Test discrete sine transform inverts", "[]") {
  Constants c("test_constants.json");

  Variable var(c, 1, true);

  for(int ix=0; ix<var.nX; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix*c.dx;
      var.spatialData[var.calcIndex(ix, k)] = 0.0;
      for(int n=0; n<c.nN; ++n) {
        var.spatialData[var.calcIndex(ix, k)] += (k+1)*sin(M_PI * n * x/c.aspectRatio);
      }
    }
  }

  var.toSpectral();
  var.toPhysical();

  for(int ix=1; ix<var.nX-1; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix*c.dx;
      real sum = 0.0;
      for(int n=0; n<c.nN; ++n) {
        sum += (k+1)*sin(M_PI * n * x/c.aspectRatio);
      }
      //cout << ix << " " << k << endl;
      require_within_error(var.spatialData[var.calcIndex(ix,k)], sum, 1e-10);
    }
  }
}

TEST_CASE("Test forward discrete sine transform", "[]") {
  Constants c("test_constants.json");

  Variable var(c, 1, true);

  for(int ix=0; ix<var.nX; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix*c.dx;
      var.spatialData[var.calcIndex(ix, k)] = 0.0;
      for(int n=0; n<c.nN; ++n) {
        var.spatialData[var.calcIndex(ix, k)] += (k+1)*sin(M_PI * n * x/c.aspectRatio);
      }
    }
  }

  var.toSpectral();

  for(int k=0; k<c.nZ; ++k) {
    for(int n=1; n<c.nN; ++n) {
      require_equal(var(n, k), k+1);
    }
  }
}

TEST_CASE("Test backward discrete sine transform", "[]") {
  Constants c("test_constants.json");

  Variable var(c, 1, true);

  for(int k=0; k<c.nZ; ++k) {
    real z = k/(c.nZ-1.0);
    for(int n=1; n<c.nN; ++n) {
      var(n,k) = n*sin(M_PI*z);
    }
  }

  var.toPhysical();

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

  Variable var(c, 1, false);

  for(int ix=0; ix<var.nX; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix*c.aspectRatio/(var.nX-1.0);
      var.spatialData[var.calcIndex(ix, k)] = cos(M_PI * (k+1) * x/c.aspectRatio) + 2.0*cos(M_PI*(k+2)*x/c.aspectRatio);
    }
  }

  var.toSpectral();
  var.toPhysical();

  for(int ix=0; ix<var.nX; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix*c.aspectRatio/(var.nX-1.0);
      require_equal(var.spatialData[var.calcIndex(ix, k)], cos(M_PI * (k+1) * x/c.aspectRatio) + 2.0*cos(M_PI*(k+2)*x/c.aspectRatio));
    }
  }
}

TEST_CASE("Test forward discrete cosine transform", "[]") {
  Constants c("test_constants.json");

  Variable var(c, 1, false);

  for(int ix=0; ix<var.nX; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix*c.aspectRatio/(var.nX-1.0);
      var.spatialData[var.calcIndex(ix, k)] = 5.0 + cos(M_PI * (k+1) * x/c.aspectRatio) + 2.0*cos(M_PI*(k+2)*x/c.aspectRatio);
    }
  }

  var.toSpectral();

  for(int k=0; k<c.nZ; ++k) {
    require_equal(var(0, k), 5.0);
    for(int n=1; n<k+1; ++n) {
      check_equal(var(n, k), 0.0);
    }

    require_equal(var((k+1), k), 1.0);
    require_equal(var((k+2), k), 2.0);

    for(int n=k+3; n<c.nN; ++n) {
      //cout << n << " " << k << endl;
      require_equal(var(n, k), 0.0);
    }
  }
}

TEST_CASE("Test backward discrete cosine transform", "[]") {
  Constants c("test_constants.json");

  Variable var(c, 1, false);

  for(int k=0; k<c.nZ; ++k) {
    var(0, k) = 5.0;
    var(2, k) = 7.0;
  }

  var.toPhysical();

  for(int ix=0; ix<var.nX; ++ix) {
    for(int k=0; k<c.nZ; ++k) {
      real x = ix/(var.nX-1.0);
      require_equal(var.spatialData[var.calcIndex(ix, k)], 5.0 + 7.0*cos(M_PI*2*x));
    }
  }
}

TEST_CASE("Test discrete Fourier transform", "[]") {
  Constants c("test_constants_periodic.json");

  Variable var(c);

  for(int k=0; k<c.nZ; ++k) {
    var(0, k) = 5.0 + k;
    var(2, k) = 7.0 + 2.0i;
    var(3, k) = 2.0 - 10.0i;
  }

  var.toPhysical();

  for(int k=0; k<c.nZ; ++k) {
    for(int ix=0; ix<var.nX; ++ix) {
      //cout << ix << " " << k << endl;
      real x = real(ix)/var.nX;
      check_equal(var.spatial(ix,k),
          5.0 + k +
          2.0*7.0*cos(2.0*M_PI*2*x) - 2.0*2.0*sin(2.0*M_PI*2*x) +
          2.0*2.0*cos(2.0*M_PI*3*x) + 2.0*10.0*sin(2.0*M_PI*3*x));
    }
  }

  var.toSpectral();

  for(int k=0; k<c.nZ; ++k) {
    check_equal(var(0, k), 5.0 + k);
    check_equal(var(2, k), 7.0 + 2.0i);
    check_equal(var(3, k), 2.0 - 10.0i);
  }
}

TEST_CASE("Test complex poisson solver", "[]") {
  Constants c("test_constants_periodic.json");

  Sim sim(c);
  Variable psi(c);

  for(int k=0; k<c.nZ; ++k) {
    for(int n=0; n<c.nN; ++n) {
      real z = c.dz*k;
      real x = c.dx*n;
      psi(n,k) = 2.0*sin(M_PI*z)*(n+1) + 3.0i*sin(2.0*M_PI*z);
    }
  }

  for(int n=0; n<c.nN; ++n) {
    sim.vars.omg(n,0) = sim.vars.omg(n,c.nZ-1) = 0.0;
    for(int k=1; k<c.nZ-1; ++k) {
      sim.vars.omg(n,k) = -(psi.dfdz2(n,k) - pow(real(n)*c.wavelength, 2)*psi(n,k));
    }
  }

  sim.solveForPsi();

  for(int k=0; k<c.nZ; ++k) {
    for(int n=0; n<c.nN; ++n) {
      //cout << n << " " << k << endl;
      require_within_error(sim.vars.psi(n,k), psi(n,k));
    }
  }
}

TEST_CASE("Test simple nonlinear multiplication and transform", "[]") {
  Constants c("test_constants.json");

  Variable var1(c, 1, false);

  for(int k=0; k<c.nZ; ++k) {
    var1(0, k) = 5.0;
    var1(2, k) = 2.0;
  }

  Variable var2(c, 1, false);

  for(int k=0; k<c.nZ; ++k) {
    var2(0, k) = 2.0;
  }

  var1.toPhysical();
  var2.toPhysical();

  Variable result(c, 1, false);

  for(int k=0; k<c.nZ; ++k) {
    for(int i=0; i<result.nX; ++i) {
      result.spatial(i,k) = 
        var1.spatial(i,k) *
        var2.spatial(i,k);
    }
  }

  result.toSpectral();

  for(int k=0; k<c.nZ; ++k) {
    require_equal(result(0,k), 2.0*var1(0,k));
    require_equal(result(2,k), 2.0*var1(2,k));
  }
}

TEST_CASE("Test spatial dfdz derivative", "[]") {
  Constants c("test_constants.json");

  Variable var(c);

  for(int i=0; i<var.nX; ++i) {
    for(int k=0; k<c.nZ; ++k) {
      real z = k/(c.nZ-1.0);
      real x = i/(var.nX-1.0);
      var.spatialData[var.calcIndex(i,k)] = sin(M_PI*z)*sin(M_PI*x);
    }
  }

  for(int k=-1; k<=c.nZ; ++k) {
    real z = k*c.dz;
    var(1,k) = sin(M_PI*z);
  }

  for(int k=1; k<c.nZ-1; ++k) {
    real z = k*c.dz;
    require_within_error(var.dfdz(1,k), M_PI*cos(M_PI*z), 1e-2);
  }

  for(int i=1; i<var.nX-1; ++i) {
    for(int k=1; k<c.nZ-1; ++k) {
      //cout << i << " " << k << endl;
      real z = k/(c.nZ-1.0);
      real x = i/(var.nX-1.0);
      require_within_error(var.dfdzSpatial(i,k), M_PI*cos(M_PI*z)*sin(M_PI*x), 1e-2);
    }
  }
}

TEST_CASE("Test spatial dfdx derivative", "[]") {
  Constants c("test_constants.json");

  Variable var(c);

  for(int i=0; i<var.nX; ++i) {
    for(int k=0; k<c.nZ; ++k) {
      real z = k/(c.nZ-1.0);
      real x = i/(var.nX-1.0);
      var.spatialData[var.calcIndex(i,k)] = sin(M_PI*z)*sin(M_PI*x);
    }
  }

  for(int i=1; i<var.nX-1; ++i) {
    for(int k=1; k<c.nZ-1; ++k) {
      real z = k/(c.nZ-1.0);
      real x = i/(var.nX-1.0);
      require_within_error(var.dfdx(i,k), M_PI/c.aspectRatio*sin(M_PI*z)*cos(M_PI*x));
    }
  }
}

TEST_CASE("Test spatial nonlinear derivative", "[]") {
  Constants c("test_constants_periodic.json");

  Sim sim(c);

  for(int i=0; i<c.nX; ++i) {
    for(int k=0; k<c.nZ; ++k) {
      real z = k*c.dz;
      real x = i*c.dx;
      sim.vars.psi.spatial(i,k) = (pow(x,2)+2.0*x-3.0)*(pow(z,2) + 5.0*z + 7.0);
      sim.vars.omg.spatial(i,k) = (pow(x,2)/3.0+x-2.0)*(2.0*pow(z,2) + 3.0*z + 10.0);
    }
  }

  Variable nonlinearSpectralTerm(c);

  sim.computeNonlinearDerivative(nonlinearSpectralTerm, sim.vars.omg);

  for(int i=1; i<c.nX-1; ++i) {
    for(int k=1; k<c.nZ-1; ++k) {
      real z = k*c.dz;
      real x = i*c.dx;
      real dpsidx = (2.0*x + 2.0)*(pow(z,2) + 5.0*z + 7.0);
      real dpsidz = (pow(x,2)+2.0*x-3.0)*(2.0*z + 5.0);
      real domgdx = (2.0/3.0*x + 1.0)*(2.0*pow(z,2) + 3.0*z + 10.0);
      real domgdz = (pow(x,2)/3.0+x-2.0)*(4.0*z + 3.0);

      REQUIRE(sim.nonlinearSineTerm.spatial(i,k)
        == Approx(-(-dpsidz*domgdx + dpsidx*domgdz)).margin(0.5));
    }
  }
}

TEST_CASE("Test reading and writing from file", "[]") {
  Constants c;
  c.nN = 5;
  c.nZ = 10;
  c.nX = c.nN*3+1;
  c.calculateDerivedConstants();

  Variable tmp(c);

  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      tmp(n,k) = n+1 + k;
    }
  }

  std::ofstream outfile ("./test.dat", std::ios::out | std::ios::binary);
  tmp.writeToFile(outfile);
  outfile.close();

  Variable tmp2(c);

  std::ifstream infile ("./test.dat", std::ios::in | std::ios::binary);
  tmp2.readFromFile(infile);
  infile.close();

  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      require_equal(tmp(n,k), tmp2(n,k));
    }
  }
}
