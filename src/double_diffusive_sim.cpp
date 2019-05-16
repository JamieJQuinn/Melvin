#include <iostream>
#include <cmath>

#include <double_diffusive_sim.hpp>
#include <utility.hpp>

using namespace std;

DoubleDiffusiveSimulation::DoubleDiffusiveSimulation(const Constants &c_in)
  : Sim(c_in)
  , xi(c_in)
  , dXidt(c_in)
{}

void DoubleDiffusiveSimulation::save() {
  std::ofstream file (c.saveFolder+std::string("vars")+strFromNumber(saveNumber++)+std::string(".dat"), std::ios::out | std::ios::binary);
  if(file.is_open()) {
    tmp.writeToFile(file);
    omg.writeToFile(file);
    psi.writeToFile(file);
    dTmpdt.writeToFile(file);
    dOmgdt.writeToFile(file);
    xi.writeToFile(file);
    dXidt.writeToFile(file);
  } else {
    cout << "Couldn't open " << c.saveFolder << " for writing. Aborting." << endl;
    exit(-1);
  }
  file.close();
}

void DoubleDiffusiveSimulation::load(const std::string &icFile) {
  std::ifstream file (c.icFile, std::ios::in | std::ios::binary);
  if(file.is_open()) {
    tmp.readFromFile(file);
    omg.readFromFile(file);
    psi.readFromFile(file);
    dTmpdt.readFromFile(file);
    dOmgdt.readFromFile(file);
    xi.readFromFile(file);
    dXidt.readFromFile(file);
  } else {
    cout << "Couldn't open " << c.icFile << " for reading. Aborting." << endl;
    exit(-1);
  }
}

void DoubleDiffusiveSimulation::reinit() {
  Sim::reinit();
  xi.fill(0.0);
  dXidt.fill(0.0);
}

void DoubleDiffusiveSimulation::computeLinearDerivatives() {
  // Computes the (linear) derivatives of xi and omg
  Sim::computeLinearDerivatives();
  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      dXidt(n,k) = c.tau*(xi.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*xi(n,k));
      dOmgdt(n,k) += -c.RaXi*c.tau*c.Pr*(n*M_PI/c.aspectRatio)*xi(n,k);
    }
  }
}

void DoubleDiffusiveSimulation::addAdvectionApproximation() {
  Sim::addAdvectionApproximation();
  // Only applies to the linear simulation
  for(int k=1; k<c.nZ-1; ++k) {
    dXidt(0,k) = 0.0;
  }
  for(int n=1; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      dXidt(n,k) += -1*c.xiGrad*n*M_PI/c.aspectRatio * psi(n,k);
    }
  }
}

void DoubleDiffusiveSimulation::initialLinearConditions() {
  Sim::initialLinearConditions();
  for(int n=1; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      xi(n,k) = 0.1*sin(M_PI * k*c.dz);
    }
  }
}

void DoubleDiffusiveSimulation::updateVars(real f) {
  tmp.update(dTmpdt, dt, f);
  omg.update(dOmgdt, dt, f);
  xi.update(dXidt, dt, f);
}

void DoubleDiffusiveSimulation::advanceDerivatives() {
  dTmpdt.advanceTimestep();
  dOmgdt.advanceTimestep();
  dXidt.advanceTimestep();
}
