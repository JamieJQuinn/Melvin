#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include <precision.hpp>
#include <constants.hpp>
#include <variable.hpp>

#include <iostream>

using namespace std;

TEST_CASE( "Test periodic boundary conditions", "[]" ) {
  Constants c;
  c.nN = 5;
  c.nZ = 10;
  c.aspectRatio = 1;
  c.boundaryConditions="periodic";
  c.calculateDerivedConstants();

  // Create GPU variables
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
