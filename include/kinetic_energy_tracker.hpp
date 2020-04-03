#pragma once

#include <vector>

#include <precision.hpp>
#include <constants.hpp>
#include <variable.hpp>

class KineticEnergyTracker {
  public:
    KineticEnergyTracker(const Constants &c_in);
    void calcKineticEnergy(Variable &psi);
    real calcKineticEnergyForMode(const Variable &psi, int n);
    void saveKineticEnergy();
    real calcKineticEnergyPhysical(Variable& psi);
    real calcKineticEnergySpectral(const Variable& psi);
    void calcKineticEnergyDensity(const Variable &psi);

  private:
    // Kinetic Energy tracker
    std::vector<real> kineticEnergies;
    Variable keDens;
    const Constants c;
};
