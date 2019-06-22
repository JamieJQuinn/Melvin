#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cassert>

#include <sim.hpp>
#include <precision.hpp>
#include <utility.hpp>
#include <numerical_methods.hpp>
#include <derivatives.hpp>

using namespace std;

Sim::Sim(const Constants &c_in)
  : c(c_in)
  , psi(c_in)
  , omg(c_in)
  , tmp(c_in)
  , dTmpdt(c_in, 2)
  , dOmgdt(c_in, 2)
  , xi(c_in)
  , dXidt(c_in, 2)
{
  saveNumber=0;

  variableList.push_back(&tmp);
  variableList.push_back(&omg);
  variableList.push_back(&psi);
  variableList.push_back(&dTmpdt);
  variableList.push_back(&dOmgdt);

  if(c.isDoubleDiffusion) {
    variableList.push_back(&xi);
    variableList.push_back(&dXidt);
  }

  for(int i=0; i<variableList.size(); ++i) {
    variableList[i]->initialiseData();
  }

  dt = c.initialDt;

  thomasAlgorithm = new ThomasAlgorithm(c.nZ, c.nN, c.aspectRatio, c.oodz2);
}

Sim::~Sim() {
  delete thomasAlgorithm;
}

std::string Sim::createSaveFilename() {
  // Format save number
  char buff[10];
  sprintf(buff, "%04d", saveNumber);

  return c.saveFolder+std::string("dump")+std::string(buff)+std::string(".dat");
}

void Sim::save() {
  std::ofstream file (createSaveFilename(), std::ios::out | std::ios::binary);
  if(file.is_open()) {
    for(int i=0; i<variableList.size(); ++i) {
      variableList[i]->writeToFile(file);
    }
  } else {
    cout << "Couldn't open " << c.saveFolder << " for writing. Aborting." << endl;
    exit(-1);
  }
  file.close();

  ++saveNumber;
}

void Sim::load(const std::string &icFile) {
  std::ifstream file (c.icFile, std::ios::in | std::ios::binary);
  if(file.is_open()) {
    for(int i=0; i<variableList.size(); ++i) {
      variableList[i]->readFromFile(file);
    }
  } else {
    cout << "Couldn't open " << c.icFile << " for reading. Aborting." << endl;
    exit(-1);
  }
}

void Sim::reinit() {
  for(int i=0; i<variableList.size(); ++i) {
    variableList[i]->fill(0.0);
  }
}

void Sim::saveKineticEnergy() {
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

real Sim::calcKineticEnergyForMode(int n) {
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

void Sim::calcKineticEnergy() {
  // Uses trapezeoid rule to calc kinetic energy for each mode
  real ke = 0.0;
  for(int n=0; n<c.nN; ++n) {
    ke += calcKineticEnergyForMode(n);
  }

  kineticEnergies.push_back(ke);
}

void Sim::solveForPsi(){
  // Solve for Psi using Thomas algorithm
  for(int n=0; n<c.nN; ++n) {
    thomasAlgorithm->solve(psi.getMode(n), omg.getMode(n), n);
  }
}

void Sim::applyBoundaryConditions() {
  // Streamfunction
  for(int n=0; n<c.nN; ++n) {
    assert(psi(n,0) == 0.0);
    assert(psi(n,c.nZ-1) == 0.0);
  }
  for(int k=0; k<c.nZ; ++k) {
    assert(psi(0,k) < EPSILON);
  }

  // Temp
  for(int n=0; n<c.nN; ++n) {
    if(n>0) {
      assert(tmp(n, 0) < EPSILON);
    } else {
      assert(tmp(n, 0) - 1.0 < EPSILON);
    }
    assert(tmp(n, c.nZ-1) < EPSILON);
  }

  // Vorticity
  for(int n=0; n<c.nN; ++n) {
    // check BCs
    assert(omg(n, 0) < EPSILON);
    assert(omg(n, c.nZ-1) < EPSILON);
  }
}

void Sim::printBenchmarkData() const {
  cout << t << " of " << c.totalTime << "(" << t/c.totalTime*100 << ")" << endl;
  for(int n=0; n<21; ++n) {
    printf("%d | %e | %e | %e\n", n, tmp(n,32), omg(n,32), psi(n,32));
  }
}

void Sim::updateVars(real f) {
  tmp.update(dTmpdt, dt, f);
  omg.update(dOmgdt, dt, f);
  if(c.isDoubleDiffusion) {
    xi.update(dTmpdt, dt, f);
  }
}

void Sim::advanceDerivatives() {
  dTmpdt.advanceTimestep();
  dOmgdt.advanceTimestep();
  if(c.isDoubleDiffusion) {
    dXidt.advanceTimestep();
  }
}

void Sim::runNonLinearStep(real f) {
  computeLinearDerivatives(dTmpdt, tmp, dOmgdt, omg, dXidt, xi, c);
  computeNonlinearDerivatives(dTmpdt, tmp, dOmgdt, omg, dXidt, xi, psi, c);
  updateVars(f);
  advanceDerivatives();
  solveForPsi();
}

void Sim::runNonLinear() {
  // Load initial conditions
  load(c.icFile);

  real saveTime = 0;
  real KEcalcTime = 0;
  real KEsaveTime = 0;
  real CFLCheckTime = 0;
  real f = 1.0f; // Fractional change in dt (if CFL condition being breached)
  t = 0;
  while (c.totalTime-t>EPSILON) {
    if(KEcalcTime-t < EPSILON) {
      calcKineticEnergy();
      KEcalcTime += 1e2*dt;
    }
    if(KEsaveTime-t < EPSILON) {
      saveKineticEnergy();
      KEsaveTime += 1e4*dt;
    }
    if(CFLCheckTime-t < EPSILON) {
      cout << "Checking CFL" << endl;
      CFLCheckTime += 1e4*dt;
      f = checkCFL(psi, c.dz, c.dx, dt, c.aspectRatio, c.nN, c.nX, c.nZ);
    }
    if(saveTime-t < EPSILON) {
      cout << t << " of " << c.totalTime << "(" << t/c.totalTime*100 << "%)" << endl;
      saveTime+=c.timeBetweenSaves;
      save();
    }
    runNonLinearStep(f);
    t+=dt;
    f=1.0f;
  } 
  printf("%e of %e (%.2f%%)\n", t, c.totalTime, t/c.totalTime*100);
  save();
  saveKineticEnergy();
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
  reinit();
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
  load(c.icFile);

  // Stuff for critical rayleigh check
  real tmpPrev = tmp(nCrit, 32);
  real psiPrev = psi(nCrit, 32);
  real omgPrev = omg(nCrit, 32);

  real logTmpPrev = 0.0, logOmgPrev = 0.0, logPsiPrev = 0.0;

  real tolerance = 1e-8;
  int steps = 0;
  t=0;
  while (t<c.totalTime) {
    if(steps%500 == 0) {
      real logTmp = std::log(std::abs(tmp(nCrit,32))) - std::log(std::abs(tmpPrev));
      real logOmg = std::log(std::abs(omg(nCrit,32))) - std::log(std::abs(omgPrev));
      real logPsi = std::log(std::abs(psi(nCrit,32))) - std::log(std::abs(psiPrev));
      if(std::abs(logTmp - logTmpPrev)<tolerance and
         std::abs(logOmg - logOmgPrev)<tolerance and
         std::abs(logPsi - logPsiPrev)<tolerance) {
        return logTmp;
      }
      logTmpPrev = logTmp;
      logOmgPrev = logOmg;
      logPsiPrev = logPsi;
      tmpPrev = tmp(nCrit, 32);
      psiPrev = psi(nCrit, 32);
      omgPrev = omg(nCrit, 32);
    }
    if(steps%1000 == 0) {
      calcKineticEnergy();
    }
    steps++;
    runLinearStep();
    t+=dt;
  }
  saveKineticEnergy();
  return 0;
}

void Sim::runLinearStep() {
  computeLinearDerivatives(dTmpdt, tmp, dOmgdt, omg, dXidt, xi, c);
  addAdvectionApproximation(dTmpdt, tmp, dOmgdt, omg, dXidt, xi, psi, c);
  updateVars();
  advanceDerivatives();
  solveForPsi();
}
