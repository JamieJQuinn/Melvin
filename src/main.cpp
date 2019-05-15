#include <iostream>
#include <string>
#include <cmath>

#include <sim.hpp>
#include <precision.hpp>

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

  Sim simulation(c);

#ifdef NONLINEAR
  cout << "NONLINEAR" << endl;
  simulation.runNonLinear();
#endif

#ifdef LINEAR
  cout << "LINEAR" << endl;
  real initialRa = simulation.c.Ra;
  real RaCrits [10];
  int nCritAnalytical = c.aspectRatio/sqrt(2);
  cout << "Critial mode should be " << nCritAnalytical << endl;
  cout << "Critial Ra should be " << pow(M_PI/c.aspectRatio, 4) * pow(pow(nCritAnalytical,2) + pow(c.aspectRatio,2), 3) / pow(nCritAnalytical,2) << endl;
  for(int n=1; n<11; ++n){
    cout << "Finding critical Ra for n=" << n << endl;
    real RaLower = 0.0;
    real RaUpper = initialRa;
    while(std::abs(RaLower - RaUpper) > 1) {
      simulation.reinit();
      simulation.c.Ra = (RaUpper+RaLower)/2;
      cout << "Trying Ra=" << simulation.c.Ra << endl;
      bool isCritical = simulation.isCritical(n);
      if(simulation.isFinished()) {
        cout << "Total time breached." << endl;
        break;
      }
      if(isCritical) {
#ifdef DDC
        // Let's raise Ra to damp xi instability
        RaLower = simulation.c.Ra;
#endif
#ifndef DDC
        // Upper bound is critical, let's lower it
        RaUpper = simulation.c.Ra;
#endif
      } else {
#ifdef DDC
        // Let's lower Ra to enhance xi instability
        RaUpdder = simulation.c.Ra;
#endif
#ifndef DDC
        // Lower bound is critical, let's raise it
        RaLower = simulation.c.Ra;
#endif
      }
    }

#ifdef DDC
    cout << "Critical Ra for n=" << n << " is RaXi - Ra =" << c.RaXi - simulation.c.Ra << endl;
#endif
#ifndef DDC
    cout << "Critical Ra for n=" << n << " is Ra=" << simulation.c.Ra << endl;
#endif
    RaCrits[n-1] = simulation.c.Ra;
  }
  for(int n=1; n<11; ++n) {
    cout << RaCrits[n-1] << "\\" << endl;
  }
#endif

  cout << "ENDING SIMULATION" << endl;
  return 0;
}

