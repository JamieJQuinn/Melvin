#include <kinetic_energy_tracker.hpp>

#include <cmath>

KineticEnergyTracker::KineticEnergyTracker(const Constants &c_in):
  c(c_in),
  keDens(c_in)
{
  keDens.initialiseData();
}

void KineticEnergyTracker::calcKineticEnergyDensity(const Variable &psi) {
  int nX = c.nX;
  for(int k=0; k<c.nZ; ++k) {
    for(int i=0; i<nX; ++i) {
      keDens.spatial(i,k) = pow(psi.dfdzSpatial(i,k), 2) + pow(psi.dfdx(i,k), 2);
    }
  }
}

real KineticEnergyTracker::calcKineticEnergyForMode(const Variable &psi, int n) {
  real z0 = 0.0; // limits of integration
  real z1 = 1.0;
  real ke = 0; // Kinetic energy

  ke += pow(n*M_PI/c.aspectRatio*psi.magnitude(n,0), 2)/2.0; // f(0)/2
  ke += pow(n*M_PI/c.aspectRatio*psi.magnitude(n,c.nZ-1), 2)/2.0; // f(1)/2
  for(int k=1; k<c.nZ-1; ++k) {
    ke += pow(std::abs(psi.dfdz(n,k)), 2) + pow(n*M_PI/c.aspectRatio*psi.magnitude(n,k), 2);
  }

  ke *= (z1-z0)*c.aspectRatio/(4*(c.nZ-1));
  return ke;
}

real KineticEnergyTracker::calcKineticEnergySpectral(const Variable& psi) {
  // Uses trapezeoid rule to calc kinetic energy for each mode
  real ke = 0.0;
  for(int n=0; n<c.nN; ++n) {
    ke += calcKineticEnergyForMode(psi, n);
  }

  return ke;
}

real KineticEnergyTracker::calcKineticEnergyPhysical(Variable& psi) {
  int nX = c.nX;
  int nZ = c.nZ;
  real ke = 0.0;

  psi.toPhysical();
  calcKineticEnergyDensity(psi);

  ke += keDens.spatial(0,0);
  ke += keDens.spatial(0,nZ-1);
  ke += keDens.spatial(nX-1,0);
  ke += keDens.spatial(nX-1,nZ-1);

  for(int k=1; k<nZ-1; ++k) {
    ke += 2.0*(keDens.spatial(0,k) + keDens.spatial(nX-1,k));
  }

  for(int i=1; i<nX-1; ++i) {
    ke += 2.0*(keDens.spatial(i,0) + keDens.spatial(i,nZ-1));
  }

  for(int k=1; k<nZ-1; ++k) {
    for(int i=1; i<nX-1; ++i) {
      ke += 4.0*keDens.spatial(i,k);
    }
  }

  ke *= c.dz*c.dx/4.0;

  return ke;
}

void KineticEnergyTracker::calcKineticEnergy(Variable &psi) {
  real ke = calcKineticEnergyPhysical(psi);

  kineticEnergies.push_back(ke);
}

void KineticEnergyTracker::saveKineticEnergy() {
  // Save total energy
  std::ofstream file (c.saveFolder+"kinetic_energy.dat", std::ios::out | std::ios::binary);
  for(int i=0; i<kineticEnergies.size(); ++i) {
    file.write(reinterpret_cast<char*>(&kineticEnergies[i]), sizeof(real));
  }
  file.flush();
  file.close();
}

