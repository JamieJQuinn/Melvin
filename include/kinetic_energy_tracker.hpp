#pragma once

#include <vector>

#include <precision.hpp>
#include <constants.hpp>
#include <variable.hpp>

class KineticEnergyTracker {
  public:
    KineticEnergyTracker(const Constants &c_in);
    void calcKineticEnergy(const Variable &psi);
    real calcKineticEnergyForMode(const Variable &psi, int n);
    void saveKineticEnergy();

  private:
    // Kinetic Energy tracker
    std::vector<real> kineticEnergies;
    const Constants c;
};
