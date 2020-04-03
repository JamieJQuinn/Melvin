#pragma once

#include <precision.hpp>
#include <variable.hpp>

class ThomasAlgorithm {
  private:
    void formTriDiagonalArraysForN(const real *sub, const mode *dia, const real *sup,
        mode * wk1, mode *wk2);
    void precalculate();

    void solveSystem(mode *sol, const mode *rhs, const int matrixN, const int n) const;
    void solveSystem(Variable& sol, const Variable& rhs, const int matrixN, const int n) const;
    void solvePeriodicSystem(Variable& sol, const Variable& rhs, const int n) const;

    int nZ;
    const int nN;
    const bool isPeriodic;
    const real oodz2;
    const mode wavelength;

    // For periodic solver
    mode *sol2;
    mode *rhs2;
  public:
    mode *wk1;
    mode *wk2;
    real *sub;

    void solve(Variable& sol, const Variable& rhs, const int n) const;
    ThomasAlgorithm(const Constants& c_in);
    ~ThomasAlgorithm();
};
