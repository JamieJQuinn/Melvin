#include <kinetic_energy_tracker.hpp>

#include <cmath>

KineticEnergyTracker::KineticEnergyTracker(const Constants &c_in):
  c(c_in)
{}

real KineticEnergyTracker::calcKineticEnergyForMode(const Variable &psi, int n) {
  real z0 = 0.0; // limits of integration
  real z1 = 1.0;
  real ke = 0; // Kinetic energy
    ke += pow(n*M_PI/c.aspectRatio*psi(n,0), 2)/2.0; // f(0)/2
    ke += pow(n*M_PI/c.aspectRatio*psi(n,c.nZ-1), 2)/2.0; // f(1)/2
    for(int k=1; k<c.nZ-1; ++k) {
      ke += pow(psi.dfdz(n,k), 2) + pow(n*M_PI/c.aspectRatio*psi(n,k), 2);
    }
  ke *= (z1-z0)*c.aspectRatio/(4*(c.nZ-1));
  return ke;
}

void KineticEnergyTracker::calcKineticEnergy(const Variable &psi) {
  // Uses trapezeoid rule to calc kinetic energy for each mode
  real ke = 0.0;
  for(int n=0; n<c.nN; ++n) {
    ke += calcKineticEnergyForMode(psi, n);
  }

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
  // save energy per mode
  //for(int n=1; n<c.nN; ++n) {
    //std::ofstream file (c.saveFolder+"KineticEnergyMode"+strFromNumber(n)+std::string(".dat"), std::ios::out | std::ios::app | std::ios::binary);
    //real ke = calcKineticEnergyForMode(n);
    //file.write(reinterpret_cast<char*>(&ke), sizeof(real));
    //file.flush();
    //file.close();
  //}
}

