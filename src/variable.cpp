#include <variable.hpp>
#include <precision.hpp>
#include <numerical_methods.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <omp.h>

using namespace std::complex_literals;

mode* Variable::getCurrent() {
  return getPlus(0);
}

const mode* Variable::getCurrent() const {
  return getPlus(0);
}

mode* Variable::getPrevious() {
  return data + previous*varSize();
}

const mode* Variable::getPrevious() const {
  return data + previous*varSize();
}

mode* Variable::getPlus(int nSteps) {
  return data + ((current+nSteps)%totalSteps)*varSize();
}

const mode* Variable::getPlus(int nSteps) const {
  return data + ((current+nSteps)%totalSteps)*varSize();
}

void Variable::advanceTimestep(int nSteps) {
  // Advance the timestep on (i.e. current++)
  previous = current;
  current = (current + nSteps)%totalSteps;
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

void Variable::fill(mode value) {
  for(int i=0; i<this->totalSize(); ++i) {
    data[i] = value;
    spatialData[i] = value.real();
  }
}

void Variable::update(const Variable& dVardt, const real dt, const real f) {
  for(int k=0; k<nZ; ++k) {
    for(int n=0; n<nN; ++n) {
      (*this)(n, k) += adamsBashforth(dVardt(n, k), dVardt.getPrev(n, k), f, dt);
    }
  }
}

void Variable::initialiseData(mode initialValue) {
  data = new mode[this->totalSize()];
  spatialData = new real[this->totalSize()];
  fill(initialValue);
}

void Variable::toSpectral() {
  fftw_execute(fftwForwardPlan);

  if(c.horizontalBoundaryConditions == BoundaryConditions::impermeable) {
    #pragma omp parallel for schedule(dynamic)
    for(int k=0; k<nZ; ++k) {
      (*this)(0,k) = ((*this)(0,k))/(2.0*(nX-1.0));
      for(int n=1; n<nX; ++n) {
        (*this)(n,k) = ((*this)(n,k))/(nX-1.0);
      }
    }
  } else if(c.horizontalBoundaryConditions == BoundaryConditions::periodic) {
    #pragma omp parallel for schedule(dynamic)
    for(int k=0; k<nZ; ++k) {
      for(int n=0; n<nX; ++n) {
        (*this)(n,k) *= 1.0/nX;
      }
    }
  }
}

void Variable::toPhysical() {
  fftw_execute(fftwBackwardPlan);

  if(c.horizontalBoundaryConditions == BoundaryConditions::impermeable) {
    if(useSinTransform) {
      #pragma omp parallel for schedule(dynamic)
      for(int k=0; k<nZ; ++k) {
        for(int i=0; i<nX; ++i) {
          spatial(i,k) = (spatial(i,k))/2.0;
        }
      }
    } else {
      #pragma omp parallel for schedule(dynamic)
      for(int k=0; k<nZ; ++k) {
        for(int i=0; i<nX; ++i) {
          spatial(i,k) = (spatial(i,k) + (*this)(0,k).real() + pow(-1, i)*(*this)(nX-1,k).real())/2.0;
        }
      }
    }
  }
}

void Variable::setupFFTW() {
  if(c.horizontalBoundaryConditions == BoundaryConditions::impermeable) {
    int n[1];
    fftw_r2r_kind kind[1];
    real *spatial;
    mode *spectral;

    if(useSinTransform) {
      kind[0] = FFTW_RODFT00;
      n[0] = nX-2;
      spatial = spatialData + calcIndex(0,0) + 1;
      spectral = getCurrent() + calcIndex(0,0) + 1;
    } else {
      kind[0] = FFTW_REDFT00;
      n[0] = nX;
      spatial = spatialData + calcIndex(0,0);
      spectral = getCurrent() + calcIndex(0,0);
    }

    fftw_plan_with_nthreads(omp_get_max_threads());

    fftwForwardPlan = fftw_plan_many_r2r(1, n, nZ,
        spatial, NULL, 1, rowSize(),
        (real*)spectral, NULL, 2, 2*rowSize(),
        kind, FFTW_MEASURE);

    fftw_plan_with_nthreads(omp_get_max_threads());

    fftwBackwardPlan = fftw_plan_many_r2r(1, n, nZ,
        (real*)spectral, NULL, 2, 2*rowSize(),
        spatial, NULL, 1, rowSize(),
        kind, FFTW_MEASURE);
  } else if(c.horizontalBoundaryConditions == BoundaryConditions::periodic) {
    int n[] = {nX};
    real *spatial = spatialData + calcIndex(0,0);
    mode *spectral = getCurrent() + calcIndex(0,0);

    fftw_plan_with_nthreads(omp_get_max_threads());

    fftwForwardPlan = fftw_plan_many_dft_r2c(1, n, nZ,
        spatial, NULL, 1, rowSize(),
        (fftw_complex*)spectral, NULL, 1, rowSize(),
        FFTW_MEASURE);

    fftw_plan_with_nthreads(omp_get_max_threads());

    fftwBackwardPlan = fftw_plan_many_dft_c2r(1, n, nZ,
        (fftw_complex*)spectral, NULL, 1, rowSize(),
        spatial, NULL, 1, rowSize(),
        FFTW_MEASURE | FFTW_PRESERVE_INPUT);
  }
}

Variable::Variable(const Constants &c_in, const int totalSteps_in, const bool useSinTransform_in):
  data(nullptr),
  spatialData(nullptr),
  nN(c_in.nN),
  nX(c_in.nX),
  nZ(c_in.nZ),
  dz(c_in.dz),
  boundaryConditions(c_in.verticalBoundaryConditions),
  oodz2(c_in.oodz2),
  oodz(c_in.oodz),
  oodx(c_in.oodx),
  totalSteps(totalSteps_in),
  current(0),
  previous(1),
  nG(c_in.nG),
  useSinTransform(useSinTransform_in),
  c(c_in)
{
  initialiseData(0.0);

  #pragma omp critical
  {
    setupFFTW();
  }

  if(c.horizontalBoundaryConditions == BoundaryConditions::periodic) {
    xDerivativeFactor = 1.0i*c.wavelength;
  } else {
    if(useSinTransform) {
      xDerivativeFactor = c.wavelength;
    } else {
      xDerivativeFactor = -c.wavelength;
    }
  }
}

Variable::~Variable() {
  if(data != nullptr) {
    delete [] data;
  }
  if(spatialData != nullptr) {
    delete [] spatialData;
  }
}
