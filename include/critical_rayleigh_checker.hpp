#pragma once

#include <constants.hpp>
#include <precision.hpp>
#include <sim.hpp>
#include <sim_gpu.hpp>

#include <cmath>

class CriticalRayleighChecker {
  public:
    CriticalRayleighChecker(const Constants &c_in);
    int calculateCriticalWavenumber();
    real calculateCriticalRayleigh();
    int testCriticalRayleigh();
    template<class SimType> bool isCritical(const real testRa, const int nCrit);

  private:
    Constants c;
    bool didTestFinish;
};

template<class SimType>
bool CriticalRayleighChecker::isCritical(const real testRa, const int nCrit) {
  c.Ra = testRa;
  SimType sim(c);

  sim.vars.load(c.icFile);

  // Stuff for critical rayleigh check
  real tmpPrev = sim.vars.tmp(nCrit, 32);
  real psiPrev = sim.vars.psi(nCrit, 32);
  real omgPrev = sim.vars.omg(nCrit, 32);

  real logTmpPrev = 0.0, logOmgPrev = 0.0, logPsiPrev = 0.0;

  real tolerance = 1e-8;
  int steps = 0;
  real t=0;
  while (t<c.totalTime) {
    if(steps%500 == 0) {
      real logTmp = std::log(std::abs(sim.vars.tmp(nCrit,32))) - std::log(std::abs(tmpPrev));
      real logOmg = std::log(std::abs(sim.vars.omg(nCrit,32))) - std::log(std::abs(omgPrev));
      real logPsi = std::log(std::abs(sim.vars.psi(nCrit,32))) - std::log(std::abs(psiPrev));
      if(std::abs(logTmp - logTmpPrev)<tolerance and
         std::abs(logOmg - logOmgPrev)<tolerance and
         std::abs(logPsi - logPsiPrev)<tolerance) {
        return logTmp > 0.0;
      }
      logTmpPrev = logTmp;
      logOmgPrev = logOmg;
      logPsiPrev = logPsi;
      tmpPrev = sim.vars.tmp(nCrit, 32);
      psiPrev = sim.vars.psi(nCrit, 32);
      omgPrev = sim.vars.omg(nCrit, 32);
    }
    if(steps%1000 == 0) {
      sim.keTracker.calcKineticEnergy(sim.vars.psi);
    }
    steps++;
    sim.runLinearStep();
    t+=sim.dt;
  }
  sim.keTracker.saveKineticEnergy();
  didTestFinish = false;
  return false;
}
