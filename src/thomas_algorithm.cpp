#include <thomas_algorithm.hpp>
#include <cassert>
#include <cmath>

ThomasAlgorithm::~ThomasAlgorithm() {
  delete[] wk1;
  delete[] wk2;
  delete[] sub;
  if(sol2 != nullptr) {
    delete [] sol2;
  }
  if(rhs2 != nullptr) {
    delete [] rhs2;
  }
}

ThomasAlgorithm::ThomasAlgorithm(const int nZ_in, const int nN_in, const int a_in, const real oodz2, bool isPeriodic_in) :
  nZ {nZ_in},
  nN {nN_in},
  a {a_in},
  isPeriodic {isPeriodic_in},
  oodz2 {oodz2},
  sol2 {nullptr},
  rhs2 {nullptr}
{
  wk1 = new real [nN*nZ];
  wk2 = new real [nN*nZ];
  sub = new real [nZ];

  if(isPeriodic) {
    sol2 = new real[nZ-1];
    rhs2 = new real[nZ-1];
    for(int k=0; k<nZ-1; ++k) {
      rhs2[k] = 0.0;
    }
  }

  precalculate();
}

void ThomasAlgorithm::precalculate() {
  // Precalculate tridiagonal arrays
  real * dia = new real [nZ];
  real * sup = new real [nZ];
  for(int k=0; k<nZ; ++k) {
    sub[k] = sup[k] = -oodz2;
  }
  for(int n=0; n<nN; ++n) {
    for(int k=0; k<nZ; ++k){
      dia[k] = pow(M_PI/a*n, 2) + 2*oodz2;
    }
    if(not isPeriodic) {
      // This encodes impermeable vertical boundary conditions
      dia[0] = dia[nZ-1] = 1.0;
      sub[nZ-2] = sup[0] = 0.0;
    }
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
          const real *sub, const real *dia, const real *sup,
    real * wk1, real *wk2) {
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

void ThomasAlgorithm::solveSystem(real *sol, const real *rhs, const int matrixN, const int n) const {
  // Solves the tridiagonal system represented by sub, dia and sup.
  // If sub, dia and sup do not change, they can be rolled into wk1 and wk2
  // using formTridiArrays() and simply saved

  int iN = n*nZ;

  // Forward Subsitution
  sol[0] = rhs[0]*wk1[0+iN];
  for (int i=1; i<matrixN; ++i) {
    sol[i] = (rhs[i] - sub[i-1]*sol[i-1])*wk1[i+iN];
  }
  // Backward Substitution
  for (int i=matrixN-2; i>=0; --i) {
    sol[i] -= wk2[i+iN]*sol[i+1];
  }
}

void ThomasAlgorithm::solveSystem(Variable& sol, const Variable& rhs, const int matrixN, const int n) const {
  // Solves the tridiagonal system represented by sub, dia and sup.
  // If sub, dia and sup do not change, they can be rolled into wk1 and wk2
  // using formTridiArrays() and simply saved

  int iN = n*nZ;

  // Forward Subsitution
  sol(n,0) = rhs(n,0)*wk1[0+iN];
  for (int i=1; i<matrixN; ++i) {
    sol(n,i) = (rhs(n,i) - sub[i-1]*sol(n,i-1))*wk1[i+iN];
  }
  // Backward Substitution
  for (int i=matrixN-2; i>=0; --i) {
    sol(n,i) -= wk2[i+iN]*sol(n,i+1);
  }
}

void ThomasAlgorithm::solvePeriodicSystem(Variable& sol, const Variable& rhs, const int n) const {
  solveSystem(sol, rhs, nZ-1, n);

  real a, b, c;
  a = c = -oodz2;
  b = pow(M_PI/a*n, 2) + 2*oodz2;
  rhs2[0] = -a;
  rhs2[nZ-2] = -c;
  solveSystem(sol2, rhs2, nZ-1, n);

  real x_last = (rhs(n,nZ-1) - c*sol(n,0) - a*sol(n,nZ-2))/(b + a*sol2[nZ-2] + c*sol2[0]);

  for(int k=0; k<nZ-1; ++k) {
    sol(n,k) += x_last*sol2[k];
  }
  sol(n,nZ-1) = x_last;
}

void ThomasAlgorithm::solve(Variable& sol, const Variable& rhs, const int n) const {
  if(isPeriodic) {
    solvePeriodicSystem(sol, rhs, n);
  } else {
    solveSystem(sol, rhs, nZ, n);
  }
}
