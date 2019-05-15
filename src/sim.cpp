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
  , psi(c_in)
  , omg(c_in)
  , tmp(c_in)
#ifdef DDC
  , xi(c_in)
  , dXidt(c_in, 2)
#endif
  , dTmpdt(c_in, 2)
  , dOmgdt(c_in, 2)
{
  kePrev = 0.0f;
  keCurrent = 0.0f;
  saveNumber=0;
  KEsaveNumber=0;

  dt = c.initialDt;

  thomasAlgorithm = new ThomasAlgorithm(c.nZ, c.nN, c.aspectRatio, c.oodz2);
}

Sim::~Sim() {
  delete thomasAlgorithm;
}

void Sim::save() {
  std::ofstream file (c.saveFolder+std::string("vars")+strFromNumber(saveNumber++)+std::string(".dat"), std::ios::out | std::ios::binary);
  if(file.is_open()) {
    tmp.writeToFile(file);
    omg.writeToFile(file);
    psi.writeToFile(file);
    dTmpdt.writeToFile(file);
    dOmgdt.writeToFile(file);
  } else {
    cout << "Couldn't open " << c.saveFolder << " for writing. Aborting." << endl;
    exit(-1);
  }
  file.close();
}

void Sim::load(const std::string &icFile) {
  std::ifstream file (c.icFile, std::ios::in | std::ios::binary);
  if(file.is_open()) {
    tmp.readFromFile(file);
    omg.readFromFile(file);
    psi.readFromFile(file);
    dTmpdt.readFromFile(file);
    dOmgdt.readFromFile(file);
  } else {
    cout << "Couldn't open " << c.icFile << " for reading. Aborting." << endl;
    exit(-1);
  }
}

void Sim::reinit() {
  tmp.fill(0.0);
  omg.fill(0.0);
  psi.fill(0.0);
  dTmpdt.fill(0.0);
  dOmgdt.fill(0.0);
#ifdef DDC
  xi.fill(0.0);
  dXidt.fill(0.0);
#endif
}

void Sim::saveKineticEnergy() {
  // Save total energy
  std::ofstream file (c.saveFolder+"KineticEnergy"+std::string(".dat"), std::ios::out | std::ios::app | std::ios::binary);
  real ke = calcKineticEnergy();
  kePrev = keCurrent;
  keCurrent = ke;
  file.write(reinterpret_cast<char*>(&ke), sizeof(real));
  file.flush();
  file.close();
  // save energy per mode
  for(int n=1; n<c.nN; ++n) {
    std::ofstream file (c.saveFolder+"KineticEnergyMode"+strFromNumber(n)+std::string(".dat"), std::ios::out | std::ios::app | std::ios::binary);
    real ke = calcKineticEnergyForMode(n);
    file.write(reinterpret_cast<char*>(&ke), sizeof(real));
    file.flush();
    file.close();
  }
}

//#ifdef DDC
//void Sim::updateXi(real f=1.0) {
  //for(int n=0; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //xi[n*c.nZ+k] += adamsBashforth(dXidt[current*c.nZ*c.nN+n*c.nZ+k], dXidt[((current+1)%2)*c.nZ*c.nN+n*c.nZ+k], f, dt);
    //}
  //} 
//}
//#endif

void Sim::updateTmpAndOmg(real f = 1.0) {
  // Update variables using Adams-Bashforth Scheme
  // f is the proportional change between the new dt and old dt
  // ( if dt changed )
  for(int n=0; n<c.nN; ++n) {
    for(int k=0; k<c.nZ; ++k) {
      tmp(n, k) += adamsBashforth(dTmpdt(n, k), dTmpdt.getPrev(n, k), f, dt);
      omg(n, k) += adamsBashforth(dOmgdt(n, k), dOmgdt.getPrev(n, k), f, dt);

      assert(!isnan(tmp(n, k)));
      assert(!isnan(omg(n, k)));
    }
    // check BCs
    if(n>0) {
      assert(tmp(n, 0) < EPSILON);
    } else {
      assert(tmp(n, 0) - 1.0 < EPSILON);
    }
    assert(tmp(n, c.nZ-1) < EPSILON);
    assert(omg(n, 0) < EPSILON);
    assert(omg(n, c.nZ-1) < EPSILON);
  }

  // Boundary Conditions
  // Periodic
  //tmp[3*c.nZ+0] = sin(OMEGA*t);
}


void Sim::computeLinearDerivatives(int linearSim) {
  // Computes the (linear) derivatives of Tmp and omg
  // If linear sim is 0, we start n from 0 and the advection approximation
  // in dTmpdt vanishes
  for(int n=linearSim; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      dTmpdt(n,k) = tmp.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*tmp(n,k);
#ifdef DDC
      dXidt(n,k) = c.tau*(xi.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*xi(n,k));
#endif
      if (linearSim == 1) {
#ifdef DDC
        dXidt(n,k) += -1*c.xiGrad*n*M_PI/c.aspectRatio * psi(n,k);
#endif
        dTmpdt(n,k) += -1*c.tempGrad*n*M_PI/c.aspectRatio * psi(n,k);
      }
      assert(!isnan(tmp.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*tmp(n,k)
        + n*M_PI/c.aspectRatio * psi(n,k)*linearSim));
      dOmgdt(n,k) =
        c.Pr*(omg.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*omg(n,k)
        + c.Ra*n*M_PI/c.aspectRatio*tmp(n,k)
        );
#ifdef DDC
      dOmgdt(n,k) += -c.c.RaXi*c.tau*c.Pr*(n*M_PI/c.aspectRatio)*xi(n,k);
#endif
      assert(dOmgdt(0,k) < EPSILON);
    }
  }
}

void Sim::computeNonLinearDerivatives() {
  for(int n=1; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      // Contribution TO tmp[n=0]
      dTmpdt(0,k) +=
        -M_PI/(2*c.aspectRatio)*n*(
          psi.dfdz(n,k)*tmp(n,k) +
          tmp.dfdz(n,k)*psi(n,k)
          );
    }
  }
  #pragma omp parallel for schedule(dynamic)
  for(int n=1; n<c.nN; ++n) {
    // Contribution FROM tmp[n=0]
    for(int k=1; k<c.nZ-1; ++k) {
      dTmpdt(n,k) +=
        -n*M_PI/c.aspectRatio*psi(n,k)*tmp.dfdz(0,k);
    }
    // Contribution FROM tmp[n>0] and omg[n>0]
    int o;
    for(int m=1; m<n; ++m){
      // Case n = n' + n''
      o = n-m;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dTmpdt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          -m*psi.dfdz(o,k)*tmp(m,k)
          +o*tmp.dfdz(m,k)*psi(o,k)
          );
        dOmgdt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          -m*psi.dfdz(o,k)*omg(m,k)
          +o*omg.dfdz(m,k)*psi(o,k)
          );
      }
    }
    for(int m=n+1; m<c.nN; ++m){
      // Case n = n' - n''
      o = m-n;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dTmpdt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          +m*psi.dfdz(o,k)*tmp(m,k)
          +o*tmp.dfdz(m,k)*psi(o,k)
          );
        dOmgdt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          +m*psi.dfdz(o,k)*omg(m,k)
          +o*omg.dfdz(m,k)*psi(o,k)
          );
      }
    }
    for(int m=1; m+n<c.nN; ++m){
      // Case n= n'' - n'
      o = n+m;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dTmpdt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          +m*psi.dfdz(o,k)*tmp(m,k)
          +o*tmp.dfdz(m,k)*psi(o,k)
          );
        dOmgdt(n,k) +=
          +M_PI/(2*c.aspectRatio)*(
          +m*psi.dfdz(o,k)*omg(m,k)
          +o*omg.dfdz(m,k)*psi(o,k)
          );
      }
    }
  }
}

void Sim::solveForPsi(){
  // Solve for Psi using Thomas algorithm
  for(int n=0; n<c.nN; ++n) {
    thomasAlgorithm->solve(psi.getMode(n), omg.getMode(n), n);
    // Check Boundary Conditions
    assert(psi(n,0) == 0.0);
    assert(psi(n,c.nZ-1) == 0.0);
  }
  // Check BCs
  for(int k=0; k<c.nZ; ++k) {
    assert(psi(0,k) < EPSILON);
  }

}

void Sim::printBenchmarkData() const {
  cout << t << " of " << c.totalTime << "(" << t/c.totalTime*100 << ")" << endl;
  for(int n=0; n<21; ++n) {
    printf("%d | %e | %e | %e\n", n, tmp(n,32), omg(n,32), psi(n,32));
  }
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

real Sim::calcKineticEnergy() {
  // Uses trapezeoid rule to calc kinetic energy for each mode
  real ke = 0.0;
  for(int n=0; n<c.nN; ++n) {
    ke += calcKineticEnergyForMode(n);
  }
  return ke;
}

void Sim::runNonLinear() {
  // Load initial conditions
  load(c.icFile);
  real saveTime = 0;
  real KEsaveTime = 0;
  real CFLCheckTime = 0;
  real f = 1.0f; // Fractional change in dt (if CFL condition being breached)
  t = 0;
  while (c.totalTime-t>EPSILON) {
    if(KEsaveTime-t < EPSILON) {
      saveKineticEnergy();
      KEsaveTime += 1e-4;
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
    computeLinearDerivatives(0);
    computeNonLinearDerivatives();
    updateTmpAndOmg(f);
    f=1.0f;
    solveForPsi();
    t+=dt;
    dTmpdt.advanceTimestep();
    dOmgdt.advanceTimestep();
  } 
  printf("%e of %e (%.2f%%)\n", t, c.totalTime, t/c.totalTime*100);
  save();
}

real Sim::findCriticalRa(int nCrit) {
  load(c.icFile);

  // Stuff for critical rayleigh check
  real tmpPrev = tmp(nCrit, 32);
  real psiPrev = psi(nCrit, 32);
  real omgPrev = omg(nCrit, 32);
#ifdef DDC
  real xiPrev = xi(nCrit, 32);
#endif

  real logTmpPrev = 0.0, logOmgPrev = 0.0, logPsiPrev = 0.0;
#ifdef DDC
  real logXiPrev = 0.0;
#endif

  real tolerance = 1e-10;
  int steps = 0;
  t=0;
  while (t<c.totalTime) {
    if(steps%500 == 0) {
      real logTmp = std::log(std::abs(tmp(nCrit,32))) - std::log(std::abs(tmpPrev));
#ifdef DDC
      real logXi = std::log(std::abs(xi(nCrit,32))) - std::log(std::abs(xiPrev));
#endif
      real logOmg = std::log(std::abs(omg(nCrit,32))) - std::log(std::abs(omgPrev));
      real logPsi = std::log(std::abs(psi(nCrit,32))) - std::log(std::abs(psiPrev));
      if(std::abs(logTmp - logTmpPrev)<tolerance and
#ifdef DDC
         std::abs(logXi  -  logXiPrev)<tolerance and
#endif
         std::abs(logOmg - logOmgPrev)<tolerance and
         std::abs(logPsi - logPsiPrev)<tolerance) {
        return logTmp;
      }
      logTmpPrev = logTmp;
#ifdef DDC
      logXiPrev = logXi;
#endif
      logOmgPrev = logOmg;
      logPsiPrev = logPsi;
      tmpPrev = tmp(nCrit, 32);
      psiPrev = psi(nCrit, 32);
      omgPrev = omg(nCrit, 32);
#ifdef DDC
      xiPrev = xi(nCrit, 32);
#endif
    }

    steps++;
    computeLinearDerivatives(1);
    updateTmpAndOmg();
#ifdef DDC
    updateXi();
#endif
    solveForPsi();
    t+=dt;
    dTmpdt.advanceTimestep();
    dOmgdt.advanceTimestep();
  }
  return 0;
}
