#include <iostream>
#include <string>
#include <cmath>

#include <sim.hpp>
#include <precision.hpp>
#include <double_diffusive_sim.hpp>

#define strVar(variable) #variable
#define OMEGA 2*M_PI*4

using namespace std;

int main(int argc, char** argv) {
  cout <<"STARTING SIMULATION\n" << endl;

  std::string constantsFile = "";
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--constants") {
      constantsFile = argv[++i];
    }
  }

  if(constantsFile == "") {
    std::cout << "No constants file found. Aborting." << std::endl;
    return -1;
  }

  const Constants c(constantsFile);
  if(not c.isValid()) {
    cout << "ABORTING\n" << endl;
    return -1;
  }
  c.print();

#ifndef DDC
  Sim simulation(c);
#endif

#ifdef DDC
  DoubleDiffusiveSimulation simulation(c);
#endif

#ifdef NONLINEAR
  cout << "NONLINEAR" << endl;
  simulation.runNonLinear();
#endif

#ifdef LINEAR
  cout << "LINEAR" << endl;
  int nCritAnalytical = c.aspectRatio/sqrt(2) + 0.5;
  real Ra_mn = pow(M_PI/c.aspectRatio, 4) * pow(pow(nCritAnalytical,2) + pow(c.aspectRatio,2), 3) / pow(nCritAnalytical,2);
  real RaCrit = Ra_mn;

#ifdef DDC
  if(c.tempGrad > 0) {
    RaCrit = c.RaXi - Ra_mn;
  } else {
    RaCrit = Ra_mn;
  }
#endif

  cout << "Critical mode should be " << nCritAnalytical << endl;
  cout << "Corresponding Ra_mn is " << Ra_mn << endl;
  cout << "And RaCrit is " << RaCrit << endl;

  simulation.c.Ra = RaCrit - 2;
  cout << "Testing Ra = " << simulation.c.Ra << endl;
  bool isBelowCritical = simulation.isCritical(nCritAnalytical);
  cout << "Below this, critical = " << isBelowCritical << endl;
  if(simulation.isFinished()) {
    cout << "Total time breached." << endl;
    return -1;
  }

  simulation.c.Ra = RaCrit + 2;
  cout << "Testing Ra = " << simulation.c.Ra << endl;
  simulation.reinit();
  bool isAboveCritical = simulation.isCritical(nCritAnalytical);
  cout << "Above this, critical = " << isAboveCritical << endl;
  if(simulation.isFinished()) {
    cout << "Total time breached." << endl;
    return -1;
  }

#ifndef DDC
  bool success = isAboveCritical and (not isBelowCritical);
#endif

#ifdef DDC
  bool success = (not isAboveCritical) and isBelowCritical;
#endif

  if(success) {
    cout << "Critical Ra FOUND." << endl;
    return 1;
  } else {
    cout << "Critical Ra NOT FOUND." << endl;
    return -1;
  }
#endif

  cout << "ENDING SIMULATION" << endl;
  return 0;
}

