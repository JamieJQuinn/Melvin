#include <variable.hpp>
#include <precision.hpp>
#include <numerical_methods.hpp>
#include <cassert>
#include <cmath>
#include <iostream>

real* Variable::getCurrent() {
  return getPlus(0);
}

const real* Variable::getCurrent() const {
  return getPlus(0);
}

real* Variable::getPrevious() {
  return data + previous*varSize();
}

const real* Variable::getPrevious() const {
  return data + previous*varSize();
}

real* Variable::getPlus(int nSteps) {
  return data + ((current+nSteps)%totalSteps)*varSize();
}

const real* Variable::getPlus(int nSteps) const {
  return data + ((current+nSteps)%totalSteps)*varSize();
}

real *Variable::getMode(int n) {
  // Exposes entire mode for the Thomas algorithm
  return getCurrent() + calcIndex(n, 0);
}

void Variable::advanceTimestep(int nSteps) {
  // Advance the timestep on (i.e. current++)
  previous = current;
  current = (current + nSteps)%totalSteps;
}

bool Variable::anyNan() const {
  for(int i=0; i<this->totalSize(); ++i) {
    if(std::isnan(data[i])) {
      return true;
    }
  }
  return false;
}

void Variable::writeToFile(std::ofstream& file) const {
  for(int i=0; i<totalSteps; ++i) {
    for(int n=0; n<nN; ++n) {
      file.write(reinterpret_cast<const char*>(getPlus(i)+calcIndex(n,0)), sizeof(data[0])*nZ);
    }
  }
}

void Variable::readFromFile(std::ifstream& file) {
  for(int i=0; i<totalSteps; ++i) {
    for(int n=0; n<nN; ++n) {
      file.read(reinterpret_cast<char*>(getPlus(i)+calcIndex(n,0)), sizeof(data[0])*nZ);
    }
  }
  topBoundary = (*this)(0,nZ-1);
  bottomBoundary = (*this)(0,0);
}

void Variable::fill(real value) {
  for(int i=0; i<this->totalSize(); ++i) {
    data[i] = value;
  }
}

void Variable::update(const Variable& dVardt, const real dt, const real f) {
  for(int n=0; n<nN; ++n) {
    for(int k=0; k<nZ; ++k) {
      (*this)(n, k) += adamsBashforth(dVardt(n, k), dVardt.getPrev(n, k), f, dt);
    }
  }
}

void Variable::initialiseData(real initialValue) {
  data = new real[this->totalSize()];
  fill(initialValue);
}

void Variable::applyBoundaryConditions() {
  if(boundaryConditions == BoundaryConditions::dirichlet) {
    for(int n=0; n<nN; ++n) {
      (*this)(n,0) = 0.0;
      (*this)(n,nZ-1) = 0.0;
    }
    (*this)(0,nZ-1) = topBoundary;
    (*this)(0,0) = bottomBoundary;
  } else if(boundaryConditions == BoundaryConditions::periodic) {
    for(int n=0; n<nN; ++n) {
      (*this)(n,-1) = (*this)(n, nZ-1);
      (*this)(n,nZ) = (*this)(n, 0);
    }
  }
}

Variable::Variable(const Constants &c_in, const int totalSteps_in):
  data(NULL),
  nN(c_in.nN),
  nZ(c_in.nZ),
  dz(c_in.dz),
  boundaryConditions(c_in.verticalBoundaryConditions),
  oodz2(c_in.oodz2),
  oodz(c_in.oodz),
  totalSteps(totalSteps_in),
  current(0),
  previous(1),
  nG(c_in.nG)
{}

Variable::~Variable() {
  if(data != NULL) {
    delete [] data;
  }
}
