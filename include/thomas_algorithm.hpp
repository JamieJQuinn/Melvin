#pragma once

#include <precision.hpp>

class ThomasAlgorithm {
  private:
    void formTriDiagonalArraysForN(const real *sub, const real *dia, const real *sup,
        real * wk1, real *wk2);

    const int nZ;
    const real oodz2;

    real *wk1;
    real *wk2;
    real *sub;
  public:
    void solve(real *sol, const real *rhs, const int n) const;
    ThomasAlgorithm(const int nZ, const int nN, const int a, const real oodz2);
    ~ThomasAlgorithm();
};
