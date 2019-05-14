#pragma once

#include <cassert>
#include <cmath>

class ThomasAlgorithm {
  private:
    void formTriDiagonalArraysForN(const double *sub, const double *dia, const double *sup,
        double * wk1, double *wk2);

    const int nZ;
    const double oodz2;

    double *wk1;
    double *wk2;
    double *sub;
  public:
    void solve(double *sol, const double *rhs, const int n) const;
    ThomasAlgorithm(const int nZ, const int nN, const int a, const double oodz2);
    ~ThomasAlgorithm();
};
