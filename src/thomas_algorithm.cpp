#include "thomas_algorithm.hpp"

ThomasAlgorithm::~ThomasAlgorithm() {
  delete[] wk1;
  delete[] wk2;
  delete[] sub;
}

ThomasAlgorithm::ThomasAlgorithm(const int nZ, const int nN, const int a, const double oodz2) :
nZ {nZ},
oodz2 {oodz2}
{
  wk1 = new double [nN*nZ];
  wk2 = new double [nN*nZ];
  sub = new double [nZ];

  // Precalculate tridiagonal arrays
  double * dia = new double [nZ];
  double * sup = new double [nZ];
  for(int k=0; k<nZ; ++k) {
    sub[k] = sup[k] = -oodz2;
  }
  for(int n=0; n<nN; ++n) {
    for(int k=0; k<nZ; ++k){
      dia[k] = pow(M_PI/a*n, 2) + 2*oodz2;
    }
    dia[0] = dia[nZ-1] = 1.0;
    sub[nZ-2] = sup[0] = 0.0;
    formTriDiagonalArraysForN(
    sub, dia, sup,
    wk1+n*nZ, wk2+n*nZ);
  }
  for(int i=0; i<nZ*nN; ++i) {
    assert(!std::isnan(wk1[i]));
    assert(!std::isnan(wk2[i]));
  }

  delete [] dia;
  delete [] sup;
}

void ThomasAlgorithm::formTriDiagonalArraysForN (
          const double *sub, const double *dia, const double *sup,
    double * wk1, double *wk2) {
  assert(dia[0] != 0.0);

  wk1[0] = 1.0/dia[0];
  wk2[0] = sup[0]*wk1[0];

  for (int i=1; i<nZ-1; ++i) {
    assert((dia[i] - sub[i-1] * wk2[i-1]) != 0.0);

    wk1[i] = 1.0/(dia[i] - sub[i-1] * wk2[i-1]);
    wk2[i] = sup[i]*wk1[i];
  }

  assert((dia[nZ-1] - sub[nZ-2]*wk2[nZ-2]) != 0.0);

  wk1[nZ-1] = 1.0/(dia[nZ-1] - sub[nZ-2]*wk2[nZ-2]);
}

void ThomasAlgorithm::solve(double *sol, const double *rhs, const int n) const {
  // Solves the tridiagonal system represented by sub, dia and sup.
  // If sub, dia and sup do not change, they can be rolled into wk1 and wk2
  // using formTridiArrays() and simply saved

  int iN = n*nZ;

  // Forward Subsitution
  assert(!std::isnan(wk1[0+iN]));
  assert(!std::isnan(rhs[0]));
  sol[0] = rhs[0]*wk1[0+iN];
  for (int i=1; i<nZ; ++i) {
    sol[i] = (rhs[i] - sub[i-1]*sol[i-1])*wk1[i+iN];

    assert(!std::isnan(rhs[i]));
    assert(!std::isnan(sub[i-1]));
    assert(!std::isnan(sol[i-1]));
    assert(!std::isnan(wk1[i+iN]));
    assert(!std::isnan(sol[i]));
  }
  // Backward Substitution
  for (int i=nZ-2; i>=0; --i) {
    sol[i] -= wk2[i+iN]*sol[i+1];
  }
}
