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

//void Sim::reinit() {
  //for(int i=0; i<c.nZ*c.nN; ++i) {
    //psi[i] = 0.0;
    //omg[i] = 0.0;
    //tmp[i] = 0.0;
//#ifdef DDC
    //xi[i] = 0.0;
//#endif
  //}
  //for(int i=0; i<c.nZ*c.nN*2; ++i) {
    //dTmpdt[i] = 0.0;
    //dOmgdt[i] = 0.0;
//#ifdef DDC
    //dXidt[i] = 0.0;
//#endif
  //}
//}

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
        dXidt(n,k) += -1*xiGrad*n*M_PI/c.aspectRatio * psi(n,k);
#endif
        dTmpdt(n,k) += -1*tmpGrad*n*M_PI/c.aspectRatio * psi(n,k);
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

//void Sim::printMaxOf(real *a, std::string name) {
  //int nStart = 0; // n level to start from
  //// Find max
  //real max = a[nStart*c.nZ];
  //int maxLoc[] = {0, nStart};
  //for(int n=nStart; n<c.nN; ++n) {
    //for(int k=0; k<c.nZ; ++k) {
      //if(a[n*c.nZ+k]>max) {
        //max = a[n*c.nZ+k];
        //maxLoc[0] = k;
        //maxLoc[1] = n;
      //}
    //}
  //}
  //// print max
  //cout << max << " @ " << "(" << maxLoc[0] << ", " << maxLoc[1] << ")" << endl;
//}

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
  current = 0;
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
      cout << std::log(std::abs(keCurrent)) - std::log(std::abs(kePrev)) << endl;
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

//real Sim::runLinear(int nCrit) {
  //// Initial Conditions
  //// Let psi = omg = dtmpdt = domgdt = 0
  //// Let tmp[n>0] = sin(PI*z)
  //// and tmp[n=0] = (1-z)/N
  //// For DDC Salt-fingering
  //real tmpGrad = 1;
//#ifdef DDC
  //real xiGrad = 1;
//#endif
  //[>
  //// For DDC SemiConvection
  //tmpGrad = -1;
//#ifdef DDC
  //xiGrad = -1;
//#endif
  //*/
//#ifndef DDC
  //tmpGrad = -1;
//#endif
  //for(int k=0; k<c.nZ; ++k) {
    //if(tmpGrad==-1){
      //tmp[c.nZ*0+k] = 1-k*c.dz;
    //} else if(tmpGrad==1) {
      //tmp[c.nZ*0+k] = k*c.dz;
    //}
//#ifdef DDC
    //if(xiGrad==-1){
      //xi[c.nZ*0+k] = 1-k*c.dz;
    //} else if(xiGrad==1) {
      //xi[c.nZ*0+k] = k*c.dz;
    //}
//#endif
    //for(int n=1; n<c.nN; ++n) {
      //tmp[c.nZ*n+k] = sin(M_PI*k*c.dz);
//#ifdef DDC
      //xi[c.nZ*n+k] = sin(M_PI*k*c.dz);
//#endif
    //} 
  //}
  //// Check BCs
  //for(int n=0; n<c.nN; ++n){
    //if(n>0) {
      //assert(tmp[n*c.nZ] < EPSILON);
    //} else {
      //assert(tmp[n*c.nZ] - 1.0 < EPSILON);
    //}
    //assert(tmp[n*c.nZ+c.nZ-1] < EPSILON);
    //assert(omg[n*c.nZ] < EPSILON);
    //assert(omg[n*c.nZ+c.nZ-1] < EPSILON);
  //}

  //// Stuff for critical rayleigh check
  //real tmpPrev[c.nN];
//#ifdef DDC
  //real xiPrev[c.nN];
//#endif
  //real omgPrev[c.nN];
  //real psiPrev[c.nN];
  //for(int n=0; n<c.nN; ++n){
    //tmpPrev[n] = tmp[32+n*c.nZ];
//#ifdef DDC
    //xiPrev[n]  = xi [32+n*c.nZ];
//#endif
    //psiPrev[n] = psi[32+n*c.nZ];
    //omgPrev[n] = omg[32+n*c.nZ];
  //}
  //real logTmpPrev = 0.0;
//#ifdef DDC
  //real logXiPrev = 0.0;
//#endif
  //real logPsiPrev =0.0;
  //real logOmgPrev =0.0;
  //real tolerance = 1e-10;
  //current = 0;
  //int steps = 0;
  //t=0;
  //while (t<c.totalTime) {
    //if(steps%500 == 0) {
      //real logTmp = std::log(std::abs(tmp[32+nCrit*c.nZ])) - std::log(std::abs(tmpPrev[nCrit]));
//#ifdef DDC
      //real logXi = std::log(std::abs(xi[32+nCrit*c.nZ])) - std::log(std::abs(xiPrev[nCrit]));
//#endif
      //real logOmg = std::log(std::abs(omg[32+nCrit*c.nZ])) - std::log(std::abs(omgPrev[nCrit]));
      //real logPsi = std::log(std::abs(psi[32+nCrit*c.nZ])) - std::log(std::abs(psiPrev[nCrit]));
      //if(std::abs(logTmp - logTmpPrev)<tolerance) {
//#ifdef DDC
      //if(std::abs(logXi - logXiPrev)<tolerance) {
//#endif
      //if(std::abs(logOmg - logOmgPrev)<tolerance) {
      //if(std::abs(logPsi - logPsiPrev)<tolerance) {
        //return logTmp;
//#ifdef DDC
      //}
//#endif
      //}}}
      //logTmpPrev = logTmp;
//#ifdef DDC
      //logXiPrev = logXi;
//#endif
      //logOmgPrev = logOmg;
      //logPsiPrev = logPsi;
      //for(int n=1; n<11; ++n){
        //tmpPrev[n] = tmp[32+n*c.nZ];
//#ifdef DDC
        //xiPrev[n] =  xi [32+n*c.nZ];
//#endif
        //psiPrev[n] = psi[32+n*c.nZ];
        //omgPrev[n] = omg[32+n*c.nZ];
      //}
    //}
    //steps++;
    //computeLinearDerivatives(1);
    //updateTmpAndOmg();
//#ifdef DDC
    //updateXi();
//#endif
    //solveForPsi();
    //t+=dt;
    //++current%=2;
  //} 
  //return 0;
//}
