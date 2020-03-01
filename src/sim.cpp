#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cassert>

#include <sim.hpp>
#include <precision.hpp>
#include <utility.hpp>
#include <numerical_methods.hpp>
#include <boundary_conditions.hpp>

using namespace std;

Sim::Sim(const Constants &c_in)
  : c(c_in)
  , vars(c_in)
  , nonlinearTerm(c_in)
  , keTracker(c_in)
{
  dt = c.initialDt;

  nonlinearTerm.initialiseData();

  thomasAlgorithm = new ThomasAlgorithm(c.nZ, c.nN, c.aspectRatio, c.oodz2, c.verticalBoundaryConditions == BoundaryConditions::periodic);
}

Sim::~Sim() {
  delete thomasAlgorithm;
}

void Sim::solveForPsi(){
  // Solve for Psi using Thomas algorithm
  #pragma omp parallel for schedule(dynamic)
  for(int n=0; n<c.nN; ++n) {
    thomasAlgorithm->solve(vars.psi, vars.omg, n);
  }
}

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
  for(int k=0; k<c.nZ; ++k) {
    for(int n=0; n<c.nN; ++n) {
      vars.dTmpdt(n,k) = vars.tmp.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*vars.tmp(n,k);
    }
  }
}

void Sim::computeLinearVorticityDerivative() {
  for(int k=0; k<c.nZ; ++k) {
    for(int n=0; n<c.nN; ++n) {
      vars.dOmgdt(n,k) =
        c.Pr*(vars.omg.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*vars.omg(n,k))
        + c.Pr*c.Ra*(n*M_PI/c.aspectRatio)*vars.tmp(n,k);
    }
  }
}

void Sim::computeLinearXiDerivative() {
  for(int k=0; k<c.nZ; ++k) {
    for(int n=0; n<c.nN; ++n) {
      vars.dXidt(n,k) = c.tau*(vars.xi.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*vars.xi(n,k));
      vars.dOmgdt(n,k) += -c.RaXi*c.tau*c.Pr*(n*M_PI/c.aspectRatio)*vars.xi(n,k);
    }
  }
}

void Sim::addAdvectionApproximation() {
  // Only applies to the linear simulation
  for(int k=0; k<c.nZ; ++k) {
    vars.dOmgdt(0,k) = 0.0;
    vars.dTmpdt(0,k) = 0.0;
  }
  for(int k=0; k<c.nZ; ++k) {
    for(int n=1; n<c.nN; ++n) {
      vars.dTmpdt(n,k) += -1*vars.tmp.dfdz(0,k)*n*M_PI/c.aspectRatio * vars.psi(n,k);
    }
  }
  if(c.isDoubleDiffusion) {
    for(int k=0; k<c.nZ; ++k) {
      vars.dXidt(0,k) = 0.0;
    }
    for(int k=0; k<c.nZ; ++k) {
      for(int n=1; n<c.nN; ++n) {
        vars.dXidt(n,k) += -1*vars.xi.dfdz(0,k)*n*M_PI/c.aspectRatio * vars.psi(n,k);
      }
    }
  }
}

void Sim::applyPhysicalBoundaryConditions() {
  real nX = vars.tmp.nX;
  #pragma omp parallel for schedule(dynamic)
  for(int k=0; k<c.nZ; ++k) {
    // Non-conducting
    vars.tmp.spatial(-1,k) = vars.tmp.spatial(1, k);
    vars.tmp.spatial(nX,k) = vars.tmp.spatial(nX-2, k);

    // Impermeable
    vars.psi.spatial(0,k) = 0.0;
    vars.psi.spatial(nX-1,k) = 0.0;

    // Stress free
    vars.psi.spatial(-1,k) = 2.0*vars.psi.spatial(0,k) - vars.psi.spatial(1,k);
    vars.psi.spatial(nX,k) = 2.0*vars.psi.spatial(nX-1,k) - vars.psi.spatial(nX-2,k);
    vars.omg.spatial(0,k) = 0.0;
    vars.omg.spatial(nX-1,k) = 0.0;
  }
}

void Sim::computeNonlinearDerivatives() {
  vars.tmp.toPhysical(false);
  vars.omg.toPhysical(true);
  vars.psi.toPhysical(true);
  applyPhysicalBoundaryConditions();
  computeNonlinearTemperatureDerivative();
  computeNonlinearVorticityDerivative();
  if(c.isDoubleDiffusion) {
    computeNonlinearXiDerivative();
  }
}

void Sim::computeNonlinearDerivative(Variable &dVardt, const Variable &var, const bool useSinTransform) {
  #pragma omp parallel for schedule(dynamic)
  for(int k=0; k<c.nZ; ++k) {
    for(int ix=0; ix<var.nX; ++ix) {
      nonlinearTerm.spatialData[nonlinearTerm.calcIndex(ix,k)] = 
        vars.psi.dfdzSpatial(ix, k)*var.dfdx(ix, k) -
        vars.psi.dfdx(ix, k)*var.dfdzSpatial(ix, k);
    }
  }

  nonlinearTerm.toSpectral(useSinTransform);

  for(int k=0; k<c.nZ; ++k) {
    for(int n=0; n<c.nN; ++n) {
      dVardt(n,k) += nonlinearTerm(n,k);
    }
  }
}

void Sim::computeNonlinearTemperatureDerivative() {
  computeNonlinearDerivative(vars.dTmpdt, vars.tmp, false);
}

void Sim::computeNonlinearVorticityDerivative() {
  computeNonlinearDerivative(vars.dOmgdt, vars.omg, true);
}

void Sim::computeNonlinearXiDerivative() {
  computeNonlinearDerivative(vars.dXidt, vars.xi, false);
}

void Sim::applyTemperatureBoundaryConditions() {
  vars.tmp.applyBoundaryConditions();
}

void Sim::applyVorticityBoundaryConditions() {
  for(int k=0; k<c.nZ; ++k) {
    vars.omg(0,k) = 0.0;
  }
  for(int n=0; n<c.nN; ++n) {
    vars.omg(n,0) = 0.0;
    vars.omg(n,c.nZ-1) = 0.0;
  }
}

void Sim::applyPsiBoundaryConditions() {
  for(int k=0; k<c.nZ; ++k) {
    vars.psi(0,k) = 0.0;
  }
  for(int n=0; n<c.nN; ++n) {
    vars.psi(n,0) = 0.0;
    vars.psi(n,c.nZ-1) = 0.0;
  }
}

void Sim::runNonLinearStep(real f) {
  computeLinearDerivatives();
  computeNonlinearDerivatives();
  if(c.verticalBoundaryConditions == BoundaryConditions::periodic) {
    for(int k=0; k<c.nZ; ++k) {
      vars.dTmpdt(0,k) = vars.tmp.topBoundary - vars.tmp.bottomBoundary;
    }
    if(c.isDoubleDiffusion) {
      for(int k=0; k<c.nZ; ++k) {
        vars.dXidt(0,k) = vars.xi.topBoundary - vars.xi.bottomBoundary;
      }
    }
  }
  vars.updateVars(dt, f);
  applyTemperatureBoundaryConditions();
  applyVorticityBoundaryConditions();
  vars.advanceDerivatives();
  solveForPsi();
  applyPsiBoundaryConditions();
}

void Sim::runNonLinear() {
  // Load initial conditions
  vars.load(c.icFile);

  applyTemperatureBoundaryConditions();
  applyVorticityBoundaryConditions();
  applyPsiBoundaryConditions();

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
      dt *= f;
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

void Sim::runLinearStep() {
  computeLinearDerivatives();
  addAdvectionApproximation();
  vars.updateVars(dt);
  vars.tmp.applyBoundaryConditions();
  vars.omg.applyBoundaryConditions();
  vars.advanceDerivatives();
  solveForPsi();
  vars.psi.applyBoundaryConditions();
}
