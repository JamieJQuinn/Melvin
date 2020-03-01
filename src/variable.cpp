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
      for(int k=0; k<nZ; ++k) {
        file.write(reinterpret_cast<const char*>(getPlus(i)+calcIndex(n,k)), sizeof(data[0]));
      }
    }
  }
}

void Variable::readFromFile(std::ifstream& file) {
  for(int i=0; i<totalSteps; ++i) {
    for(int n=0; n<nN; ++n) {
      for(int k=0; k<nZ; ++k) {
        file.read(reinterpret_cast<char*>(getPlus(i)+calcIndex(n,k)), sizeof(data[0]));
      }
    }
  }
  topBoundary = (*this)(0,nZ-1);
  bottomBoundary = (*this)(0,0);
}

void Variable::fill(real value) {
  for(int i=0; i<this->totalSize(); ++i) {
    data[i] = value;
    spatialData[i] = value;
  }
}

void Variable::update(const Variable& dVardt, const real dt, const real f) {
  for(int k=0; k<nZ; ++k) {
    for(int n=0; n<nN; ++n) {
      (*this)(n, k) += adamsBashforth(dVardt(n, k), dVardt.getPrev(n, k), f, dt);
    }
  }
}

void Variable::initialiseData(real initialValue) {
  data = new real[this->totalSize()];
  spatialData = new real[this->totalSize()];
  fill(initialValue);
}

void Variable::applyBoundaryConditions() {
  if(boundaryConditions == BoundaryConditions::dirichlet) {
    for(int n=0; n<nN; ++n) {
      (*this)(n,0) = 0.0;
      (*this)(n,nZ-1) = 0.0;
      // Ensure derivaties are zero
      (*this)(n,-1) = (*this)(n,1);
      (*this)(n,nZ) = (*this)(n,nZ-2);
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

void Variable::toSpectral(const bool useSinTransform) {
  fftw_plan plan;
  real *input, *output;
  if(useSinTransform) {
    plan = fftwSinePlan;
    input = spatialData + 1;
    output = getCurrent() + 1;
  } else {
    plan = fftwCosinePlan;
    input = spatialData;
    output = getCurrent();
  }

  #pragma omp parallel for schedule(dynamic)
  for(int k=0; k<nZ; ++k) {
    fftw_execute_r2r(plan, input + calcIndex(0,k), output + calcIndex(0,k));
  }
  for(int k=0; k<nZ; ++k) {
    (*this)(0,k) = ((*this)(0,k))/(2.0*(nX-1.0));
    for(int n=1; n<nX; ++n) {
      (*this)(n,k) = ((*this)(n,k))/(nX-1.0);
    }
  }
}

void Variable::toPhysical(const bool useSinTransform) {
  fftw_plan plan;
  real *input, *output;
  if(useSinTransform) {
    plan = fftwSinePlan;
    input = getCurrent() + 1;
    output = spatialData + 1;
  } else {
    plan = fftwCosinePlan;
    input = getCurrent();
    output = spatialData;
  }

  #pragma omp parallel for schedule(dynamic)
  for(int k=0; k<nZ; ++k) {
    fftw_execute_r2r(plan, input + calcIndex(0,k), output + calcIndex(0,k));
  }
  if(useSinTransform) {
    #pragma omp parallel for schedule(dynamic)
    for(int k=0; k<nZ; ++k) {
      for(int i=0; i<nX; ++i) {
        spatialData[calcIndex(i,k)] = (spatialData[calcIndex(i,k)])/2.0;
      }
    }
  } else {
    #pragma omp parallel for schedule(dynamic)
    for(int k=0; k<nZ; ++k) {
      for(int i=0; i<nX; ++i) {
        spatialData[calcIndex(i,k)] = (spatialData[calcIndex(i,k)] + (*this)(0,k) + pow(-1, i)*(*this)(nX-1,k))/2.0;
      }
    }
  }
}

Variable::Variable(const Constants &c_in, const int totalSteps_in):
  data(nullptr),
  spatialData(nullptr),
  nN(c_in.nN),
  nX(c_in.nN*3 + 1),
  nZ(c_in.nZ),
  dz(c_in.dz),
  boundaryConditions(c_in.verticalBoundaryConditions),
  oodz2(c_in.oodz2),
  oodz(c_in.oodz),
  totalSteps(totalSteps_in),
  current(0),
  previous(1),
  nG(c_in.nG)
{
  oodx = (nX-1.0)/c_in.aspectRatio;
  fftwSinePlan = fftw_plan_r2r_1d(nX-2, spatialData+1, getCurrent()+1, FFTW_RODFT00, FFTW_ESTIMATE);
  fftwCosinePlan = fftw_plan_r2r_1d(nX, spatialData, getCurrent(), FFTW_REDFT00, FFTW_ESTIMATE);
}

Variable::~Variable() {
  if(data != nullptr) {
    delete [] data;
  }
  if(spatialData != nullptr) {
    delete [] spatialData;
  }

  fftw_destroy_plan(fftwSinePlan);
  fftw_destroy_plan(fftwCosinePlan);
}
