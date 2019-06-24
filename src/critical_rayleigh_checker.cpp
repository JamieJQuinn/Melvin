#include <critical_rayleigh_checker.hpp>
#include <sim.hpp>

#include <cmath>
#include <iostream>

using namespace std;

CriticalRayleighChecker::CriticalRayleighChecker(const Constants &c_in):
  c(c_in)
{}

int CriticalRayleighChecker::calculateCriticalWavenumber() {
  return c.aspectRatio/sqrt(2) + 0.5;
}

real CriticalRayleighChecker::calculateCriticalRayleigh() {
  int nCritAnalytical = calculateCriticalWavenumber();
  real Ra_mn = pow(M_PI/c.aspectRatio, 4) * pow(pow(nCritAnalytical,2) + pow(c.aspectRatio,2), 3) / pow(nCritAnalytical,2);

  if(c.isDoubleDiffusion) {
    // This is ONLY for salt-fingering. Semiconvection not implemented
    return c.RaXi - Ra_mn;
  } else {
    return Ra_mn;
  }
}

int CriticalRayleighChecker::testCriticalRayleigh() {
  real RaCrit = calculateCriticalRayleigh();
  int nCritAnalytical = calculateCriticalWavenumber();

  cout << "Critical mode should be " << nCritAnalytical << endl;
  cout << "And RaCrit is " << RaCrit << endl;

  // Test below calculated critical number

  real testRa = RaCrit - 2;
  cout << "Testing Ra = " << testRa << endl;
  bool isBelowCritical = false;
  if(c.isCudaEnabled) {
#ifdef CUDA
    isBelowCritical = isCritical<SimGPU>(testRa, nCritAnalytical);
#endif
  } else {
    isBelowCritical = isCritical<Sim>(testRa, nCritAnalytical);
  }
  cout << "Below this, critical = " << isBelowCritical << endl;
  if(didTestFinish) {
    cout << "Total time breached." << endl;
  }

  // Test above calculated critical number

  testRa = RaCrit + 2;
  cout << "Testing Ra = " << testRa << endl;
  bool isAboveCritical = false;
  if(c.isCudaEnabled) {
#ifdef CUDA
    isAboveCritical = isCritical<SimGPU>(testRa, nCritAnalytical);
#endif
  } else {
    isAboveCritical = isCritical<Sim>(testRa, nCritAnalytical);
  }
  cout << "Above this, critical = " << isAboveCritical << endl;
  if(didTestFinish) {
    cout << "Total time breached." << endl;
  }

  bool success = false;
  if(c.isDoubleDiffusion) {
    success = (not isAboveCritical) and isBelowCritical;
  } else {
    success = isAboveCritical and (not isBelowCritical);
  }

  if(success) {
    cout << "Critical Ra FOUND." << endl;
    return 1;
  } else {
    cout << "Critical Ra NOT FOUND." << endl;
    return -1;
  }
}
