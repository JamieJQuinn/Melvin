#pragma once

#include <precision.hpp>
#include <variable.hpp>

class ThomasAlgorithm {
  private:
    void formTriDiagonalArraysForN(const real *sub, const real *dia, const real *sup,
        real * wk1, real *wk2);
    void precalculate();

    void solveSystem(real *sol, const real *rhs, const int matrixN, const int n) const;
    void solveSystem(Variable& sol, const Variable& rhs, const int matrixN, const int n) const;
    void solvePeriodicSystem(Variable& sol, const Variable& rhs, const int n) const;

    int nZ;
    const int nN;
    const int a;
    const bool isPeriodic;
    const real oodz2;

    // For periodic solver
    real *sol2;
    real *rhs2;
  public:
    real *wk1;
    real *wk2;
    real *sub;

    void solve(Variable& sol, const Variable& rhs, const int n) const;
    ThomasAlgorithm(const int nZ, const int nN, const int a, const real oodz2, bool isPeriodic);
    ~ThomasAlgorithm();
};
