#pragma once

#include <precision.hpp>
#include <constants.hpp>
#include <fstream>
#include <cassert>

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

    inline real& operator()(int n, int k);
    inline const real& operator()(int n, int k) const;

    inline real& getPrev(int n, int k);
    inline const real& getPrev(int n, int k) const;

    inline real dfdz(int n, int k) const;
    inline real dfdz2(int n, int k) const;

    void update(const Variable& dVardt, const real dt, const real f=1.0);

    bool anyNan() const;

    void fill(real value);

    void writeToFile(std::ofstream& file) const;
    void readFromFile(std::ifstream& file);

    void initialiseData(real initialValue = 0.0);

    // totalSteps gives the number of arrays to store, including the current one
    Variable(const Constants &c_in, const int totalSteps_in = 1);
    ~Variable();

    const int nN;
    const int nZ;
    real * data;
  protected:
    const double dz;
    const double oodz2;
    const int totalSteps;
    int current; // index pointing to slice of array representing current time
    int previous;
};

inline real& Variable::operator()(int n, int k) {
  // Get at n, k, possibly at previous step
  assert(n>=0);
  assert(n<nN);
  assert(k>=0);
  assert(k<nZ);
  return data[current*nN*nZ + n*nZ + k];
}

inline const real& Variable::operator()(int n, int k) const {
  // Get at n, k, possibly at previous step
  assert(n>=0);
  assert(n<nN);
  assert(k>=0);
  assert(k<nZ);
  return data[current*nN*nZ + n*nZ + k];
}

inline real& Variable::getPrev(int n, int k) {
  // Get at n, k, possibly at previous step
  assert(n>=0);
  assert(n<nN);
  assert(k>=0);
  assert(k<nZ);
  return data[previous*nN*nZ + n*nZ + k];
}

inline const real& Variable::getPrev(int n, int k) const {
  // Get at n, k, possibly at previous step
  assert(n>=0);
  assert(n<nN);
  assert(k>=0);
  assert(k<nZ);
  return data[previous*nN*nZ + n*nZ + k];
}

inline int Variable::totalSize() const {
  return nN*nZ*totalSteps;
}

inline int Variable::varSize() const {
  return nN*nZ;
}

inline real Variable::dfdz(int n, int k) const {
  // Avoid derivatives at the edge
  assert(k>0);
  assert(k<nZ-1);
  return ((*this)(n, k+1) - (*this)(n, k-1))/(2*this->dz);
}

inline real Variable::dfdz2(int n, int k) const {
  // Avoid derivatives at the edge
  assert(k>0);
  assert(k<nZ-1);
  return ((*this)(n, k+1) - 2*(*this)(n, k) + (*this)(n, k-1))*this->oodz2;
}
