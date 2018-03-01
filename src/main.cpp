#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include "sim.hpp"
#include "precision.hpp"

#define strVar(variable) #variable
#define OMEGA 2*M_PI*4

using namespace std;

int main(int argc, char** argv) {
// Sim::Sim(int nZ, int nN, double dt, double Ra, double Pr, int a ,double timeBetweenSaves, bool modifydt, int current, double t, double totalTime
  int nZ =-1;
  int nN =-1;
  int a = -1;
  double dt = -1.;
  double Ra = -1.;
#ifdef DDC
  double RaXi = -1.;
  double tau = -1.;
#endif 
  double Pr = -1.;
  double totalTime = -1;
  double saveTime = -1.;
  std::string saveFolder = "";
  std::string icFile = "";
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-nZ") {
      nZ = atoi(argv[++i]);
    } else if (arg == "-nN") {
      nN = atoi(argv[++i]);
    } else if (arg == "-a") {
      a = atoi(argv[++i]);
    } else if (arg == "-dt") {
      dt = atof(argv[++i]);
    } else if (arg == "-Ra") {
      Ra = atof(argv[++i]);
    } else if (arg == "-Pr") {
      Pr = atof(argv[++i]);
    } else if (arg == "-T") {
      totalTime = atof(argv[++i]);
    } else if (arg == "-S") {
      saveTime = atof(argv[++i]);
    } else if (arg == "-o") {
      saveFolder = argv[++i];
    } else if (arg == "-i") {
      icFile = argv[++i];
    }
#ifdef DDC
    else if (arg == "-RaXi") {
      RaXi = atof(argv[++i]);
    } else if (arg == "-tau") {
      tau =atof( argv[++i]);
    }
#endif
  }

  if(nZ <=0 or nN <=0 or a <= 0) {
    cout << " nZ (" << nZ
    << ") nN (" << nN
    << ") a (" << a
    << ") should be positive integers. Aborting.\n" << endl;
    return -1;
  }
  if(dt <= 0.0f
  or Ra <= 0.0f
#ifdef DDC
  or RaXi <= 0.0f
  or tau <= 0.0f
#endif
  or Pr <= 0.0f
  or totalTime <= 0.0f
  or saveTime <= 0.0f) {
    cout << " dt (" << dt
    << ") Ra (" << Ra
#ifdef DDC
    << ") RaXi (" << Ra
    << ") tau (" << Ra
#endif
    << ") Pr (" << Pr
    << ") totalTime (" << totalTime
    << ") saveTime (" << saveTime
    << ") should be positive decimals. Aborting.\n" << endl;
    return -1;
  }
  if(saveFolder == "" or icFile == "") {
    cout <<"Save folder and initial conditions file should be present. Aborting.\n" << endl;
    return -1;
  }

  cout <<"STARTING SIMULATION\n" << endl;

  cout <<"Parameters:" << endl;
  cout << "nZ: " << nZ << endl;
  cout << "nN: " << nN << endl;
  cout << "a: " << a << endl;
  cout << "Ra: " << Ra << endl;
#ifdef DDC
  cout << "RaXi: " << RaXi << endl;
  cout << "tau: " << tau << endl;
#endif
  cout << "Pr: " << Pr << endl;
  cout << "dt: " << dt << endl;
  cout << "totalTime: " << totalTime << endl;
  cout << "saveFolder: " << saveFolder << endl;
  cout << "icFile: " << icFile << endl;
  Sim simulation = Sim(nZ, nN, dt, Ra, Pr, a,
#ifdef DDC
           RaXi, tau,
#endif
      saveTime, false, 0, 0, totalTime, saveFolder, icFile);

#ifdef NONLINEAR
  cout << "NONLINEAR" << endl;
  simulation.runNonLinear();
#endif
#ifdef LINEAR
  cout << "LINEAR" << endl;
  double initialRa = simulation.Ra;
  double RaCrits [10];
  for(int n=1; n<11; ++n){
    cout << "Finding critical Ra for n=" << n << endl;
    double RaLower = 0.0;
    double RaUpper = initialRa;
    while(std::abs(RaLower - RaUpper) > 1e-3) {
      simulation.reinit();
      simulation.Ra = (RaUpper+RaLower)/2;
      cout << "Trying Ra=" << simulation.Ra << endl;
      double result = simulation.runLinear(n);
#ifdef DDC
      if(result > 0.0) {
        RaLower = simulation.Ra;
      } else if(result < 0.0) {
        RaUpper = simulation.Ra;
      } else {
        cout << "Total time breached." << endl;
        break;
      }
#endif
#ifndef DDC
      if(result < 0.0) {
        RaLower = simulation.Ra;
      } else if(result > 0.0) {
        RaUpper = simulation.Ra;
      } else {
        cout << "Total time breached." << endl;
        break;
      }
#endif
    }
    cout << "Critical Ra for n=" << n << " is Ra=" << simulation.Ra << endl;
    RaCrits[n-1] = simulation.Ra;
  }
  for(int n=1; n<11; ++n) {
    cout << RaCrits[n-1] << "\\" << endl;
  }
#endif

  cout << "ENDING SIMULATION" << endl;
  return 0;
}

