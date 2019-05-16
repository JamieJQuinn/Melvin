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
      dXidt(n,k) += -1*xi.dfdz(0,k)*n*M_PI/c.aspectRatio * psi(n,k);
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

void DoubleDiffusiveSimulation::runLinearStep() {
  computeLinearDerivatives();
  addAdvectionApproximation();
  updateVars();
  advanceDerivatives();
  solveForPsi();
  t+=dt;
}

void DoubleDiffusiveSimulation::runNonLinearStep(real f) {
  computeLinearDerivatives();
  computeNonLinearDerivatives();
  updateVars(f);
  advanceDerivatives();
  solveForPsi();
}

void DoubleDiffusiveSimulation::computeNonLinearDerivatives() {
  Sim::computeNonLinearDerivatives();
  for(int n=1; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      // Contribution TO xi[n=0]
      dXidt(0,k) +=
        -M_PI/(2*c.aspectRatio)*n*(
          psi.dfdz(n,k)*xi(n,k) +
          xi.dfdz(n,k)*psi(n,k)
          );
    }
  }
  #pragma omp parallel for schedule(dynamic)
  for(int n=1; n<c.nN; ++n) {
    // Contribution FROM tmp[n=0]
    for(int k=1; k<c.nZ-1; ++k) {
      dXidt(n,k) +=
        -n*M_PI/c.aspectRatio*psi(n,k)*xi.dfdz(0,k);
    }
    // Contribution FROM tmp[n>0] and omg[n>0]
    int o;
    for(int m=1; m<n; ++m){
      // Case n = n' + n''
      o = n-m;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dXidt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          -m*psi.dfdz(o,k)*xi(m,k)
          +o*xi.dfdz(m,k)*psi(o,k)
          );
      }
    }
    for(int m=n+1; m<c.nN; ++m){
      // Case n = n' - n''
      o = m-n;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dXidt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          +m*psi.dfdz(o,k)*xi(m,k)
          +o*xi.dfdz(m,k)*psi(o,k)
          );
      }
    }
    for(int m=1; m+n<c.nN; ++m){
      // Case n= n'' - n'
      o = n+m;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dXidt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          +m*psi.dfdz(o,k)*xi(m,k)
          +o*xi.dfdz(m,k)*psi(o,k)
          );
      }
    }
  }
}
