#include <variable.hpp>
#include <precision.hpp>
#include <numerical_methods.hpp>
#include <cassert>

//real* Variable::get() const {
  //return getPlus(0);
//}

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

//int Variable::writeToFile(FILE* fp, int nSteps) const {
  //int totalLength = len()*nSteps;

  //// This makes sure data is written with current step first
  //size_t nItemsWritten = 0;
  //for(int i=0; i<nSteps; ++i) {
    //nItemsWritten += std::fwrite(getPlus(i), sizeof(real), len(), fp);
  //}

  //if( int(nItemsWritten) != totalLength ) {
    //return -1;
  //} else {
    //return 0;
  //}
//}

//int Variable::writeToFile(FILE* fp) const {
  //return writeToFile(fp, totalSteps);
//}
  
//int Variable::readFromFile(FILE* fp) {
  //return readFromFile(fp, totalSteps);
//}

//void Variable::print() const {
  //printTo(std::cout);
//}

//void Variable::printTo(std::ostream& stream) const {
  //for(int i = 0; i < len(); ++i) {
    //stream << (*this)[i] << ", ";
  //}
  //stream << std::endl;
//}

//real Variable::operator[](const int i) const {
  //return data[current * len() + i];
//}

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
  for(int i=0; i<this->totalSize(); ++i) {
    data[i] = initialValue;
  }
}

Variable::~Variable() {
  delete [] data;
}
