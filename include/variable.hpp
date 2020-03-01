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
    real * getMode(int n);
    real * getCurrent();
    const real * getCurrent() const;
    real* getPrevious();
    const real* getPrevious() const;
    real * getPlus(int nSteps=0);
    const real* getPlus(int nSteps) const;

    void advanceTimestep(int nSteps=1);
    inline int totalSize() const;
    inline int varSize() const;

    inline int calcIndex(int step, int n, int k) const;
    inline int calcIndex(int n, int k) const;
    inline real& operator()(int n, int k);
    inline const real& operator()(int n, int k) const;

    inline real& getPrev(int n, int k);
    inline const real& getPrev(int n, int k) const;

    inline real dfdz(int n, int k) const;
    inline real dfdzSpatial(int ix, int k) const;
    inline real dfdx(int ix, int k) const;
    inline real dfdz2(int n, int k) const;
    inline real& spatial(int n, int k);
    inline const real& spatial(int n, int k) const;

    void update(const Variable& dVardt, const real dt, const real f=1.0);

    void applyBoundaryConditions();

    void toSpectral(const bool useSinTransform);
    void toPhysical(const bool useSinTransform);

    bool anyNan() const;

    void fill(real value);

    void writeToFile(std::ofstream& file) const;
    void readFromFile(std::ifstream& file);

    void initialiseData(real initialValue = 0.0);
    void initialiseFFTW(const bool useSinTransform);

    // totalSteps gives the number of arrays to store, including the current one
    Variable(const Constants &c_in, const int totalSteps_in = 1);
    ~Variable();

    const int nN;
    const int nZ;
    const int nX;
    const int nG;
    real * data;
    real * spatialData;

    real topBoundary;
    real bottomBoundary;
  protected:
    const BoundaryConditions boundaryConditions;
    const real dz;
    const real oodz2;
    const real oodz;
    real oodx;
    const int totalSteps;
    int current; // index pointing to slice of array representing current time
    int previous;

    fftw_plan fftwSinePlan;
    fftw_plan fftwCosinePlan;
};

inline int Variable::calcIndex(int step, int n, int k) const {
  return step*varSize() + calcIndex(n,k);
}

inline int Variable::calcIndex(int n, int k) const {
  assert(n>=0);
  assert(n<nX);
  assert(k>=-nG);
  assert(k<nZ+nG);
  return (k+nG)*(nX+2*nG) + n+nG;
  //return n*(nZ+2*nG) + (k+nG);
}

inline real& Variable::operator()(int n, int k) {
  // Get at n, k, possibly at previous step
  return data[calcIndex(current, n, k)];
}

inline const real& Variable::operator()(int n, int k) const {
  // Get at n, k, possibly at previous step
  return data[calcIndex(current, n, k)];
}

inline real& Variable::spatial(int i, int k) {
  // Get at n, k, possibly at previous step
  return spatialData[calcIndex(current, i, k)];
}

inline const real& Variable::spatial(int i, int k) const {
  // Get at n, k, possibly at previous step
  return spatialData[calcIndex(current, i, k)];
}

inline real& Variable::getPrev(int n, int k) {
  // Get at n, k, possibly at previous step
  return data[calcIndex(previous, n, k)];
}

inline const real& Variable::getPrev(int n, int k) const {
  // Get at n, k, possibly at previous step
  return data[calcIndex(previous, n, k)];
}

inline int Variable::totalSize() const {
  return varSize()*totalSteps;
}

inline int Variable::varSize() const {
  return (nX+2*nG)*(nZ+2*nG);
}

inline real Variable::dfdz(int n, int k) const {
  // Avoid derivatives at the edge
  assert(k>=0);
  assert(k<nZ);
  return ((*this)(n, k+1) - (*this)(n, k-1))*this->oodz*0.5;
}

inline real Variable::dfdzSpatial(int ix, int k) const {
  // Avoid derivatives at the edge
  assert(k>=0);
  assert(k<nZ);
  return (spatialData[calcIndex(ix, k+1)] - spatialData[calcIndex(ix, k-1)])*this->oodz*0.5;
}

inline real Variable::dfdx(int ix, int k) const {
  // Avoid derivatives at the edge
  assert(ix>0);
  assert(ix<nX-1);
  return (spatialData[calcIndex(ix+1, k)] - spatialData[calcIndex(ix-1, k)])*this->oodx*0.5;
}

inline real Variable::dfdz2(int n, int k) const {
  // Avoid derivatives at the edge
  assert(k>=0);
  assert(k<nZ);
  return ((*this)(n, k+1) - 2*(*this)(n, k) + (*this)(n, k-1))*this->oodz2;
}
