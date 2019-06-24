#include <iostream>

#include <variables.hpp>

using namespace std;

Variables::Variables(const Constants &c_in)
  : c(c_in)
  , psi(c_in)
  , omg(c_in)
  , tmp(c_in)
  , dTmpdt(c_in, 2)
  , dOmgdt(c_in, 2)
  , xi(c_in)
  , dXidt(c_in, 2)
{
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

  saveNumber = 0;
}

void Variables::updateVars(const real dt, const real f) {
  tmp.update(dTmpdt, dt, f);
  omg.update(dOmgdt, dt, f);
  if(c.isDoubleDiffusion) {
    xi.update(dTmpdt, dt, f);
  }
}

void Variables::advanceDerivatives() {
  dTmpdt.advanceTimestep();
  dOmgdt.advanceTimestep();
  if(c.isDoubleDiffusion) {
    dXidt.advanceTimestep();
  }
}

void Variables::reinit(const real value) {
  for(int i=0; i<variableList.size(); ++i) {
    variableList[i]->fill(value);
  }
}

std::string Variables::createSaveFilename() {
  // Format save number
  char buff[10];
  sprintf(buff, "%04d", saveNumber);

  return c.saveFolder+std::string("dump")+std::string(buff)+std::string(".dat");
}

void Variables::save() {
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

void Variables::load(const std::string &icFile) {
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
