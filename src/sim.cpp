#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cassert>

#include <sim.hpp>
#include <precision.hpp>
#include <utility.hpp>
#include <numerical_methods.hpp>

using namespace std;

Sim::Sim(const Constants &c_in)
  : c(c_in)
  , vars(c_in)
  , keTracker(c_in)
{
  dt = c.initialDt;

  initialiseThomasAlgorithm();
}

Sim::~Sim() {
  delete thomasAlgorithm;
}

void Sim::initialiseThomasAlgorithm() {
  thomasAlgorithm = new ThomasAlgorithm(c.nZ, c.nN, c.aspectRatio, c.oodz2);
}

void Sim::solveForPsi(){
  // Solve for Psi using Thomas algorithm
  for(int n=0; n<c.nN; ++n) {
    thomasAlgorithm->solve(vars.psi.getMode(n), vars.omg.getMode(n), n);
  }
}

////void Sim::applyBoundaryConditions() {
  //// Streamfunction
  //for(int n=0; n<c.nN; ++n) {
    //assert(vars.psi(n,0) == 0.0);
    //assert(vars.psi(n,c.nZ-1) == 0.0);
  //}
  //for(int k=0; k<c.nZ; ++k) {
    //assert(vars.psi(0,k) < EPSILON);
  //}

  //// Temp
  //for(int n=0; n<c.nN; ++n) {
    //if(n>0) {
      //assert(vars.tmp(n, 0) < EPSILON);
    //} else {
      //assert(vars.tmp(n, 0) - 1.0 < EPSILON);
    //}
    //assert(vars.tmp(n, c.nZ-1) < EPSILON);
  //}

  //// Vorticity
  //for(int n=0; n<c.nN; ++n) {
    //// check BCs
    //assert(vars.omg(n, 0) < EPSILON);
    //assert(vars.omg(n, c.nZ-1) < EPSILON);
  //}
//}

void Sim::printBenchmarkData() const {
  cout << t << " of " << c.totalTime << "(" << t/c.totalTime*100 << ")" << endl;
  for(int n=0; n<21; ++n) {
    printf("%d | %e | %e | %e\n", n, vars.tmp(n,32), vars.omg(n,32), vars.psi(n,32));
  }
}

void Sim::computeLinearDerivatives() {
  // Computes the (linear) derivatives of Tmp and vars.omg
  computeLinearTemperatureDerivative();
  computeLinearVorticityDerivative();
  if(c.isDoubleDiffusion) {
    computeLinearXiDerivative();
  }
}

void Sim::computeLinearTemperatureDerivative() {
  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      vars.dTmpdt(n,k) = vars.tmp.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*vars.tmp(n,k);
    }
  }
}

void Sim::computeLinearVorticityDerivative() {
  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      vars.dOmgdt(n,k) =
        c.Pr*(vars.omg.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*vars.omg(n,k))
        + c.Pr*c.Ra*(n*M_PI/c.aspectRatio)*vars.tmp(n,k);
    }
  }
}

void Sim::computeLinearXiDerivative() {
  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      vars.dXidt(n,k) = c.tau*(vars.xi.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*vars.xi(n,k));
      vars.dOmgdt(n,k) += -c.RaXi*c.tau*c.Pr*(n*M_PI/c.aspectRatio)*vars.xi(n,k);
    }
  }
}

void Sim::addAdvectionApproximation() {
  // Only applies to the linear simulation
  for(int k=1; k<c.nZ-1; ++k) {
    vars.dOmgdt(0,k) = 0.0;
    vars.dTmpdt(0,k) = 0.0;
  }
  for(int n=1; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      vars.dTmpdt(n,k) += -1*vars.tmp.dfdz(0,k)*n*M_PI/c.aspectRatio * vars.psi(n,k);
    }
  }
  if(c.isDoubleDiffusion) {
    for(int k=1; k<c.nZ-1; ++k) {
      vars.dXidt(0,k) = 0.0;
    }
    for(int n=1; n<c.nN; ++n) {
      for(int k=1; k<c.nZ-1; ++k) {
        vars.dXidt(n,k) += -1*vars.xi.dfdz(0,k)*n*M_PI/c.aspectRatio * vars.psi(n,k);
      }
    }
  }
}

void Sim::computeNonlinearDerivatives() {
  computeNonlinearTemperatureDerivative();
  computeNonlinearVorticityDerivative();
  if(c.isDoubleDiffusion) {
    computeNonlinearXiDerivative();
  }
}

void Sim::computeNonlinearDerivative(Variable &dVardt, const Variable &var) {
  for(int n=1; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      // Contribution TO var[n=0]
      dVardt(0,k) +=
        -M_PI/(2*c.aspectRatio)*n*(
          vars.psi.dfdz(n,k)*var(n,k) +
          var.dfdz(n,k)*vars.psi(n,k)
          );
    }
  }
  #pragma omp parallel for schedule(dynamic)
  for(int n=1; n<c.nN; ++n) {
    // Contribution FROM var[n=0]
    for(int k=1; k<c.nZ-1; ++k) {
      dVardt(n,k) +=
        -n*M_PI/c.aspectRatio*vars.psi(n,k)*var.dfdz(0,k);
    }
    // Contribution FROM var[n>0] and vars.omg[n>0]
    int o;
    for(int m=1; m<n; ++m){
      // Case n = n' + n''
      o = n-m;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dVardt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          -m*vars.psi.dfdz(o,k)*var(m,k)
          +o*var.dfdz(m,k)*vars.psi(o,k)
          );
      }
    }
    for(int m=n+1; m<c.nN; ++m){
      // Case n = n' - n''
      o = m-n;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dVardt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          +m*vars.psi.dfdz(o,k)*var(m,k)
          +o*var.dfdz(m,k)*vars.psi(o,k)
          );
      }
    }
    for(int m=1; m+n<c.nN; ++m){
      // Case n= n'' - n'
      o = n+m;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dVardt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          +m*vars.psi.dfdz(o,k)*var(m,k)
          +o*var.dfdz(m,k)*vars.psi(o,k)
          );
      }
    }
  }
}

void Sim::computeNonlinearTemperatureDerivative() {
  computeNonlinearDerivative(vars.dTmpdt, vars.tmp);
}

void Sim::computeNonlinearXiDerivative() {
  // Turns out it's the same for temperature and vars.xi
  computeNonlinearDerivative(vars.dXidt, vars.xi);
}

void Sim::computeNonlinearVorticityDerivative() {
  #pragma omp parallel for schedule(dynamic)
  for(int n=1; n<c.nN; ++n) {
    int o;
    for(int m=1; m<n; ++m){
      // Case n = n' + n''
      o = n-m;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        vars.dOmgdt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          -m*vars.psi.dfdz(o,k)*vars.omg(m,k)
          +o*vars.omg.dfdz(m,k)*vars.psi(o,k)
          );
      }
    }
    for(int m=n+1; m<c.nN; ++m){
      // Case n = n' - n''
      o = m-n;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        vars.dOmgdt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          +m*vars.psi.dfdz(o,k)*vars.omg(m,k)
          +o*vars.omg.dfdz(m,k)*vars.psi(o,k)
          );
      }
    }
    for(int m=1; m+n<c.nN; ++m){
      // Case n= n'' - n'
      o = n+m;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        vars.dOmgdt(n,k) +=
          +M_PI/(2*c.aspectRatio)*(
          +m*vars.psi.dfdz(o,k)*vars.omg(m,k)
          +o*vars.omg.dfdz(m,k)*vars.psi(o,k)
          );
      }
    }
  }
}

void Sim::runNonLinearStep(real f) {
  computeLinearDerivatives();
  computeNonlinearDerivatives();
  vars.updateVars(dt, f);
  vars.advanceDerivatives();
  solveForPsi();
}

void Sim::runNonLinear() {
  // Load initial conditions
  vars.load(c.icFile);

  real saveTime = 0;
  real KEcalcTime = 0;
  real KEsaveTime = 0;
  real CFLCheckTime = 0;
  real f = 1.0f; // Fractional change in dt (if CFL condition being breached)
  t = 0;
  while (c.totalTime-t>EPSILON) {
    if(KEcalcTime-t < EPSILON) {
      keTracker.calcKineticEnergy(vars.psi);
      KEcalcTime += 1e2*dt;
    }
    if(KEsaveTime-t < EPSILON) {
      keTracker.saveKineticEnergy();
      KEsaveTime += 1e4*dt;
    }
    if(CFLCheckTime-t < EPSILON) {
      cout << "Checking CFL" << endl;
      CFLCheckTime += 1e4*dt;
      f = checkCFL(vars.psi, c.dz, c.dx, dt, c.aspectRatio, c.nN, c.nX, c.nZ);
    }
    if(saveTime-t < EPSILON) {
      cout << t << " of " << c.totalTime << "(" << t/c.totalTime*100 << "%)" << endl;
      saveTime+=c.timeBetweenSaves;
      vars.save();
    }
    runNonLinearStep(f);
    t+=dt;
    f=1.0f;
  } 
  printf("%e of %e (%.2f%%)\n", t, c.totalTime, t/c.totalTime*100);
  vars.save();
  keTracker.saveKineticEnergy();
}

bool Sim::isCritical(int nCrit) {
  return findCriticalRa(nCrit) > 0.0;
}

real Sim::isFinished() {
  return t > c.totalTime;
}

int Sim::testCriticalRayleigh() {
  int nCritAnalytical = c.aspectRatio/sqrt(2) + 0.5;
  real Ra_mn = pow(M_PI/c.aspectRatio, 4) * pow(pow(nCritAnalytical,2) + pow(c.aspectRatio,2), 3) / pow(nCritAnalytical,2);
  real RaCrit = Ra_mn;

  if(c.isDoubleDiffusion) {
    // This is ONLY for salt-fingering. Semiconvection not implemented
    RaCrit = c.RaXi - Ra_mn;
  }

  cout << "Critical mode should be " << nCritAnalytical << endl;
  cout << "Corresponding Ra_mn is " << Ra_mn << endl;
  cout << "And RaCrit is " << RaCrit << endl;

  c.Ra = RaCrit - 2;
  cout << "Testing Ra = " << c.Ra << endl;
  bool isBelowCritical = isCritical(nCritAnalytical);
  cout << "Below this, critical = " << isBelowCritical << endl;
  if(isFinished()) {
    cout << "Total time breached." << endl;
  }

  c.Ra = RaCrit + 2;
  cout << "Testing Ra = " << c.Ra << endl;
  vars.reinit();
  bool isAboveCritical = isCritical(nCritAnalytical);
  cout << "Above this, critical = " << isAboveCritical << endl;
  if(isFinished()) {
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

real Sim::findCriticalRa(int nCrit) {
  vars.load(c.icFile);

  // Stuff for critical rayleigh check
  real tmpPrev = vars.tmp(nCrit, 32);
  real psiPrev = vars.psi(nCrit, 32);
  real omgPrev = vars.omg(nCrit, 32);

  real logTmpPrev = 0.0, logOmgPrev = 0.0, logPsiPrev = 0.0;

  real tolerance = 1e-8;
  int steps = 0;
  t=0;
  while (t<c.totalTime) {
    if(steps%500 == 0) {
      real logTmp = std::log(std::abs(vars.tmp(nCrit,32))) - std::log(std::abs(tmpPrev));
      real logOmg = std::log(std::abs(vars.omg(nCrit,32))) - std::log(std::abs(omgPrev));
      real logPsi = std::log(std::abs(vars.psi(nCrit,32))) - std::log(std::abs(psiPrev));
      if(std::abs(logTmp - logTmpPrev)<tolerance and
         std::abs(logOmg - logOmgPrev)<tolerance and
         std::abs(logPsi - logPsiPrev)<tolerance) {
        return logTmp;
      }
      logTmpPrev = logTmp;
      logOmgPrev = logOmg;
      logPsiPrev = logPsi;
      tmpPrev = vars.tmp(nCrit, 32);
      psiPrev = vars.psi(nCrit, 32);
      omgPrev = vars.omg(nCrit, 32);
    }
    if(steps%1000 == 0) {
      keTracker.calcKineticEnergy(vars.psi);
    }
    steps++;
    runLinearStep();
    t+=dt;
  }
  keTracker.saveKineticEnergy();
  return 0;
}

void Sim::runLinearStep() {
  computeLinearDerivatives();
  addAdvectionApproximation();
  vars.updateVars(dt);
  vars.advanceDerivatives();
  solveForPsi();
}
