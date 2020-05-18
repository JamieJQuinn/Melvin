#pragma once

#include <precision.hpp>
#include <constants.hpp>
#include <fstream>
#include <cassert>
#include <string>
#include <boundary_conditions.hpp>

#include <fftw3.h>

class Variable {
  // Encapsulates an array representing a variable in the model
  public:
    mode * getMode(int n);
    mode * getCurrent();
    const mode * getCurrent() const;
    mode* getPrevious();
    const mode* getPrevious() const;
    mode * getPlus(int nSteps=0);
    const mode* getPlus(int nSteps) const;

    void advanceTimestep(int nSteps=1);
    inline int totalSize() const;
    inline int varSize() const;
    inline int rowSize() const;

    inline int calcIndex(int step, int n, int k) const;
    inline int calcIndex(int n, int k) const;

    inline mode& operator()(int n, int k);
    inline const mode& operator()(int n, int k) const;
    inline real& spatial(int ix, int k);
    inline const real& spatial(int ix, int k) const;

    inline real magnitude(int n, int k) const;

    inline mode& getPrev(int n, int k);
    inline const mode& getPrev(int n, int k) const;

    inline mode dfdz(int n, int k) const;
    inline mode dfdz2(int n, int k) const;

    inline real dfdzSpatial(int ix, int k) const;
    inline real dfdx(int ix, int k) const;
    inline mode dfdxSpectral(int n, int k) const;
    inline mode dfdx2Spectral(int n, int k) const;
    inline mode laplacian(int n, int k) const;

    void update(const Variable& dVardt, const real dt, const real f=1.0);

    void toSpectral();
    void toPhysical();

    bool anyNan() const;

    void fill(mode value);

    void writeToFile(std::ofstream& file) const;
    void readFromFile(std::ifstream& file);

    void initialiseData(mode initialValue = 0.0);
    void setupFFTW();

    // totalSteps gives the number of arrays to store, including the current one
    Variable(const Constants &c_in, const int totalSteps_in = 1, const bool useSinTransform_in = true);
    ~Variable();

    const int nN;
    const int nZ;
    const int nX;
    const int nG;
    mode * data;
    real * spatialData;

    const bool useSinTransform;

    mode topBoundary;
    mode bottomBoundary;
  protected:
    const BoundaryConditions boundaryConditions;
    const real dz;
    const real oodz2;
    const real oodz;
    const real oodx;
    mode xDerivativeFactor;
    const int totalSteps;
    const Constants c;
    int current; // index pointing to slice of array representing current time
    int previous;

    fftw_plan fftwForwardPlan;
    fftw_plan fftwBackwardPlan;
};

inline int Variable::calcIndex(int step, int n, int k) const {
  return step*varSize() + calcIndex(n,k);
}

inline int Variable::calcIndex(int n, int k) const {
  return (k+nG)*(nX+2*nG) + n+nG;
  //return n*(nZ+2*nG) + (k+nG);
}

inline mode& Variable::operator()(int n, int k) {
  // Get at n, k, possibly at previous step
  return data[calcIndex(current, n, k)];
}

inline const mode& Variable::operator()(int n, int k) const {
  // Get at n, k, possibly at previous step
  return data[calcIndex(current, n, k)];
}

inline real Variable::magnitude(int n, int k) const {
  return std::abs((*this)(n,k));
}

inline real& Variable::spatial(int i, int k) {
  // Get at n, k, possibly at previous step
  return spatialData[calcIndex(current, i, k)];
}

inline const real& Variable::spatial(int i, int k) const {
  // Get at n, k, possibly at previous step
  return spatialData[calcIndex(current, i, k)];
}

inline mode& Variable::getPrev(int n, int k) {
  // Get at n, k, possibly at previous step
  return data[calcIndex(previous, n, k)];
}

inline const mode& Variable::getPrev(int n, int k) const {
  // Get at n, k, possibly at previous step
  return data[calcIndex(previous, n, k)];
}

inline int Variable::totalSize() const {
  return varSize()*totalSteps;
}

inline int Variable::varSize() const {
  return (nX+2*nG)*(nZ+2*nG);
}

inline int Variable::rowSize() const {
  return nX + 2*nG;
}

inline mode Variable::dfdz(int n, int k) const {
  // Avoid derivatives at the edge
  return ((*this)(n, k+1) - (*this)(n, k-1))*oodz*0.5;
}

inline mode Variable::dfdz2(int n, int k) const {
  // Avoid derivatives at the edge
  return ((*this)(n, k+1) - 2.0*(*this)(n, k) + (*this)(n, k-1))*oodz2;
}

inline real Variable::dfdzSpatial(int ix, int k) const {
  // Avoid derivatives at the edge
  return (spatial(ix, k+1) - spatial(ix, k-1))*oodz*0.5;
}

inline real Variable::dfdx(int ix, int k) const {
  // Avoid derivatives at the edge
  return (spatial(ix+1, k) - spatial(ix-1, k))*oodx*0.5;
}

inline mode Variable::dfdxSpectral(int n, int k) const {
  return xDerivativeFactor*real(n)*(*this)(n,k);
}

inline mode Variable::dfdx2Spectral(int n, int k) const {
  return -pow(real(n)*c.wavelength, 2)*(*this)(n,k);
}

inline mode Variable::laplacian(int n, int k) const {
  return this->dfdz2(n, k) + this->dfdx2Spectral(n, k);
}
