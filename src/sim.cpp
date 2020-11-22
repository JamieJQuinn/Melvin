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

using std::cout;
using std::endl;

Sim::Sim(const Constants &c_in)
  : c(c_in)
  , vars(c_in)
  , nonlinearSineTerm(c_in, 1, true)
  , nonlinearCosineTerm(c_in, 1, false)
  , keTracker(c_in)
{
  dt = c.initialDt;

  thomasAlgorithm = new ThomasAlgorithm(c);
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
    printf("%d | %e | %e | %e\n", n, vars.tmp.magnitude(n,32), vars.omg.magnitude(n,32), vars.psi.magnitude(n,32));
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
      vars.dTmpdt(n,k) = vars.tmp.laplacian(n,k);
    }
  }
}

void Sim::computeLinearVorticityDerivative() {
  for(int k=0; k<c.nZ; ++k) {
    for(int n=0; n<c.nN; ++n) {
      vars.dOmgdt(n,k) = c.Pr*vars.omg.laplacian(n,k) - n*c.wavelength*c.xCosDerivativeFactor*c.Pr*c.Ra*vars.tmp(n,k);
    }
  }
}

void Sim::computeLinearXiDerivative() {
  for(int k=0; k<c.nZ; ++k) {
    for(int n=0; n<c.nN; ++n) {
      vars.dXidt(n,k) = c.tau*vars.xi.laplacian(n,k);
      vars.dOmgdt(n,k) += n*c.wavelength*c.xCosDerivativeFactor*c.RaXi*c.tau*c.Pr*vars.xi(n,k);
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
      vars.dTmpdt(n,k) += -c.xSinDerivativeFactor*(n*c.wavelength)*vars.tmp.dfdz(0,k)*vars.psi(n,k);
    }
  }
  if(c.isDoubleDiffusion) {
    for(int k=0; k<c.nZ; ++k) {
      vars.dXidt(0,k) = 0.0;
    }
    for(int k=0; k<c.nZ; ++k) {
      for(int n=1; n<c.nN; ++n) {
        vars.dXidt(n,k) += -c.xSinDerivativeFactor*(n*c.wavelength)*vars.xi.dfdz(0,k)*vars.psi(n,k);
      }
    }
  }
}

void Sim::applyPhysicalBoundaryConditions() {
  real nX = c.nX;
  if(c.horizontalBoundaryConditions == BoundaryConditions::impermeable) {
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

      if(c.isDoubleDiffusion) {
        // No flux of salt
        vars.xi.spatial(-1,k) = vars.xi.spatial(1, k);
        vars.xi.spatial(nX,k) = vars.xi.spatial(nX-2, k);
      }
    }
  } else if(c.horizontalBoundaryConditions == BoundaryConditions::periodic) {
    for(int k=0; k<c.nZ; ++k) {
      vars.tmp.spatial(-1,k) = vars.tmp.spatial(nX-1, k);
      vars.tmp.spatial(nX,k) = vars.tmp.spatial(0, k);

      vars.omg.spatial(-1,k) = vars.omg.spatial(nX-1, k);
      vars.omg.spatial(nX,k) = vars.omg.spatial(0, k);

      vars.psi.spatial(-1,k) = vars.psi.spatial(nX-1, k);
      vars.psi.spatial(nX,k) = vars.psi.spatial(0, k);

      if(c.isDoubleDiffusion) {
        vars.xi.spatial(-1,k) = vars.xi.spatial(nX-1, k);
        vars.xi.spatial(nX,k) = vars.xi.spatial(0, k);
      }
    }
  }
  if (c.verticalBoundaryConditions == BoundaryConditions::periodic) {
    for(int i=0; i<c.nX; ++i) {
      vars.tmp.spatial(i, -1) = vars.tmp.spatial(i, c.nZ-1);
      vars.tmp.spatial(i, c.nZ) = vars.tmp.spatial(i, 0);

      vars.omg.spatial(i, -1) = vars.omg.spatial(i, c.nZ-1);
      vars.omg.spatial(i, c.nZ) = vars.omg.spatial(i, 0);

      vars.psi.spatial(i, -1) = vars.psi.spatial(i, c.nZ-1);
      vars.psi.spatial(i, c.nZ) = vars.psi.spatial(i, 0);

      if(c.isDoubleDiffusion) {
        vars.xi.spatial(i, -1) = vars.xi.spatial(i, c.nZ-1);
        vars.xi.spatial(i, c.nZ) = vars.xi.spatial(i, 0);
      }
    }
  }
}

void Sim::computeNonlinearDerivatives() {
  vars.tmp.toPhysical();
  vars.omg.toPhysical();
  vars.psi.toPhysical();
  if(c.isDoubleDiffusion) {
    vars.xi.toPhysical();
  }
  applyPhysicalBoundaryConditions();
  computeNonlinearTemperatureDerivative();
  computeNonlinearVorticityDerivative();
  if(c.isDoubleDiffusion) {
    computeNonlinearXiDerivative();
  }
}

void Sim::computeNonlinearDerivative(Variable &dVardt, const Variable &var) {
  Variable *nonlinearTerm;
  if(var.useSinTransform) {
    nonlinearTerm = &nonlinearSineTerm;
  } else {
    nonlinearTerm = &nonlinearCosineTerm;
  }

  #pragma omp parallel for schedule(dynamic)
  for(int k=0; k<c.nZ; ++k) {
    for(int ix=0; ix<var.nX; ++ix) {
      nonlinearTerm->spatial(ix,k) = 
        -(
            (
             var.spatial(ix+1,k)*(-vars.psi.dfdzSpatial(ix+1,k)) -
             var.spatial(ix-1,k)*(-vars.psi.dfdzSpatial(ix-1,k))
            )*c.oodx*0.5 +
            (
             var.spatial(ix,k+1)*vars.psi.dfdx(ix,k+1) -
             var.spatial(ix,k-1)*vars.psi.dfdx(ix,k-1)
            )*c.oodz*0.5
         );

    }
  }

  nonlinearTerm->toSpectral();

  for(int k=0; k<c.nZ; ++k) {
    for(int n=0; n<c.nN; ++n) {
      dVardt(n,k) += (*nonlinearTerm)(n,k);
    }
  }
}

void Sim::computeNonlinearTemperatureDerivative() {
  computeNonlinearDerivative(vars.dTmpdt, vars.tmp);
  if(c.verticalBoundaryConditions == BoundaryConditions::periodic) {
    for(int k=0; k<c.nZ; ++k) {
      vars.dTmpdt(0,k) = c.temperatureGradient;
    }
  }
}

void Sim::computeNonlinearVorticityDerivative() {
  computeNonlinearDerivative(vars.dOmgdt, vars.omg);
}

void Sim::computeNonlinearXiDerivative() {
  computeNonlinearDerivative(vars.dXidt, vars.xi);
  if(c.verticalBoundaryConditions == BoundaryConditions::periodic) {
    for(int k=0; k<c.nZ; ++k) {
      vars.dXidt(0,k) = c.salinityGradient;
    }
  }
}

void Sim::applyTemperatureBoundaryConditions() {
  if(c.verticalBoundaryConditions == BoundaryConditions::dirichlet) {
    vars.tmp(0,0) = vars.tmp.bottomBoundary;
    vars.tmp(0,c.nZ-1) = vars.tmp.topBoundary;
    for(int n=1; n<c.nN; ++n) {
      vars.tmp(n,0) = 0.0;
      vars.tmp(n,c.nZ-1) = 0.0;
    }
  } else if(c.verticalBoundaryConditions == BoundaryConditions::periodic) {
    for(int n=0; n<c.nN; ++n) {
      vars.tmp(n,-1) = vars.tmp(n,c.nZ-1);
      vars.tmp(n,c.nZ) = vars.tmp(n,0);
    }
  }
}

void Sim::applyXiBoundaryConditions() {
  if(c.verticalBoundaryConditions == BoundaryConditions::dirichlet) {
    vars.xi(0,0) = vars.xi.bottomBoundary;
    vars.xi(0,c.nZ-1) = vars.xi.topBoundary;
    for(int n=1; n<c.nN; ++n) {
      vars.xi(n,0) = 0.0;
      vars.xi(n,c.nZ-1) = 0.0;
    }
  } else if(c.verticalBoundaryConditions == BoundaryConditions::periodic) {
    for(int n=0; n<c.nN; ++n) {
      vars.xi(n,-1) = vars.xi(n,c.nZ-1);
      vars.xi(n,c.nZ) = vars.xi(n,0);
    }
  }
}

void Sim::applyVorticityBoundaryConditions() {
  if(c.verticalBoundaryConditions == BoundaryConditions::dirichlet) {
    for(int n=0; n<c.nN; ++n) {
      // Results from stress-free
      vars.omg(n,0) = 0.0;
      vars.omg(n,c.nZ-1) = 0.0;
    }
  } else if(c.verticalBoundaryConditions == BoundaryConditions::periodic) {
    for(int n=0; n<c.nN; ++n) {
      vars.omg(n,-1) = vars.omg(n,c.nZ-1);
      vars.omg(n,c.nZ) = vars.omg(n,0);
    }
  }
}

void Sim::applyPsiBoundaryConditions() {
  if(c.verticalBoundaryConditions == BoundaryConditions::dirichlet) {
    for(int n=0; n<c.nN; ++n) {
      // v_z = 0 (impermeable)
      vars.psi(n,0) = 0.0;
      vars.psi(n,c.nZ-1) = 0.0;

      // d(v_x)/dz = 0 (stress-free)
      vars.psi(n,-1) = 2.0*vars.psi(n,0) - vars.psi(n,1);
      vars.psi(n,c.nZ) = 2.0*vars.psi(n,c.nZ-1) - vars.psi(n,c.nZ-2);
    }
  } else if(c.verticalBoundaryConditions == BoundaryConditions::periodic) {
    for(int n=0; n<c.nN; ++n) {
      vars.psi(n,-1) = vars.psi(n,c.nZ-1);
      vars.psi(n,c.nZ) = vars.psi(n,0);
    }
  }
}

void Sim::runNonLinearStep(real f) {
  computeLinearDerivatives();
  computeNonlinearDerivatives();
  vars.updateVars(dt, f);
  applyTemperatureBoundaryConditions();
  applyVorticityBoundaryConditions();
  if(c.isDoubleDiffusion) {
    applyXiBoundaryConditions();
  }
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
  if(c.isDoubleDiffusion) {
    applyXiBoundaryConditions();
  }

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
      //cout << "Checking CFL" << endl;
      CFLCheckTime += 1e1*dt;
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
  applyTemperatureBoundaryConditions();
  applyVorticityBoundaryConditions();
  vars.advanceDerivatives();
  solveForPsi();
  applyPsiBoundaryConditions();
}
