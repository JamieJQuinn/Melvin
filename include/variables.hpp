#pragma once

#include <vector>
#include <iostream>

#include <variable.hpp>

template<class varType>
class Variables {
  public:
    varType tmp;
    varType xi;
    varType omg;
    varType psi;
    varType dTmpdt;
    varType dOmgdt;
    varType dXidt;

    std::vector<varType*> variableList; // for operations like saving, loading and initialising

    int saveNumber;

    const Constants c;

    Variables(const Constants &c_in);

    void save();
    void load(const std::string &icFile);
    void reinit(const real value = 0.0);

    void updateVars(const real dt, const real f = 1.0);
    void advanceDerivatives();
  private:
    std::string createSaveFilename();
};

template<class varType>
Variables<varType>::Variables(const Constants &c_in)
  : c(c_in)
  , psi(c_in, 1, true)
  , omg(c_in, 1, true)
  , tmp(c_in, 1, false)
  , dTmpdt(c_in, 2)
  , dOmgdt(c_in, 2)
  , xi(c_in, 1, false)
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

  saveNumber = 0;
}

template<class varType>
void Variables<varType>::updateVars(const real dt, const real f) {
  tmp.update(dTmpdt, dt, f);
  omg.update(dOmgdt, dt, f);
  if(c.isDoubleDiffusion) {
    xi.update(dXidt, dt, f);
  }
}

template<class varType>
void Variables<varType>::advanceDerivatives() {
  dTmpdt.advanceTimestep();
  dOmgdt.advanceTimestep();
  if(c.isDoubleDiffusion) {
    dXidt.advanceTimestep();
  }
}

template<class varType>
void Variables<varType>::reinit(const real value) {
  for(int i=0; i<variableList.size(); ++i) {
    variableList[i]->fill(value);
  }
}

template<class varType>
std::string Variables<varType>::createSaveFilename() {
  // Format save number
  char buff[10];
  sprintf(buff, "%04d", saveNumber);

  return c.saveFolder+std::string("dump")+std::string(buff)+std::string(".dat");
}

template<class varType>
void Variables<varType>::save() {
  std::ofstream file (createSaveFilename(), std::ios::out | std::ios::binary);
  if(file.is_open()) {
    for(int i=0; i<variableList.size(); ++i) {
      variableList[i]->writeToFile(file);
    }
  } else {
    std::cout << "Couldn't open " << c.saveFolder << " for writing. Aborting." << std::endl;
    exit(-1);
  }
  file.close();

  ++saveNumber;
}

template<class varType>
void Variables<varType>::load(const std::string &icFile) {
  std::ifstream file (icFile, std::ios::in | std::ios::binary);
  if(file.is_open()) {
    for(int i=0; i<variableList.size(); ++i) {
      variableList[i]->readFromFile(file);
    }
  } else {
    std::cout << "Couldn't open " << icFile << " for reading. Aborting." << std::endl;
    exit(-1);
  }
}
