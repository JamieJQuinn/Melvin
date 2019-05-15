#include <variable.hpp>
#include <precision.hpp>
#include <numerical_methods.hpp>
#include <cassert>

real* Variable::getPlus(int nSteps) {
  return data + ((current+nSteps)%totalSteps)*varSize();
}

const real* Variable::getPlus(int nSteps) const {
  return data + ((current+nSteps)%totalSteps)*varSize();
}

real *Variable::getMode(int n) {
  // Exposes entire mode for the Thomas algorithm
  return data + current*nN*nZ + n*nZ;
}

void Variable::advanceTimestep(int nSteps) {
  // Advance the timestep on (i.e. current++)
  previous = current;
  current = (current + nSteps)%totalSteps;
}


void Variable::writeToFile(std::ofstream& file) const {
  for(int i=0; i<totalSteps; ++i) {
    file.write(reinterpret_cast<const char*>(getPlus(i)), sizeof(data[0])*varSize());
  }
}

void Variable::readFromFile(std::ifstream& file) {
  for(int i=0; i<totalSteps; ++i) {
    file.read(reinterpret_cast<char*>(getPlus(i)), sizeof(data[0])*varSize());
  }
}

void Variable::fill(real value) {
  for(int i=0; i<this->totalSize(); ++i) {
    data[i] = value;
  }
}

Variable::Variable(const Constants &c_in, int totalSteps_in, real initialValue):
  data(new real[c_in.nN*c_in.nZ*totalSteps_in]()),
  nN(c_in.nN),
  nZ(c_in.nZ),
  dz(c_in.dz),
  oodz2(c_in.oodz2),
  totalSteps(totalSteps_in),
  current(0),
  previous(1)
{
  fill(initialValue);
}

Variable::~Variable() {
  delete [] data;
}
