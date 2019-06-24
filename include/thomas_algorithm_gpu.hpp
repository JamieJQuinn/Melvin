#pragma once

#include <precision.hpp>

class ThomasAlgorithmGPU {
  private:
    void formTriDiagonalArraysForN(const real *sub, const real *dia, const real *sup,
        real * wk1, real *wk2);

    const int nZ;
    const int nN;
    const real oodz2;

  public:
    real *wk1;
    real *wk2;
    real *sub;

    void solve(real *sol, const real *rhs) const;
    ThomasAlgorithmGPU(const int nZ, const int nN, const int a, const real oodz2);
    ~ThomasAlgorithmGPU();
};
