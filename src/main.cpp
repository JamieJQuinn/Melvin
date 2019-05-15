#include <iostream>
#include <string>

#include "sim.hpp"
#include "precision.hpp"

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

