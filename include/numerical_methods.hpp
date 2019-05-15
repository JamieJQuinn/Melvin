#pragma once

#include "precision.hpp"

real dfdz(real *f, int k, real dz);
real dfdz2(real *f, int k, real dz);
real adamsBashforth(real dfdt_current, real dfdt_prev, real frac, real dt);
real checkCFL(real* psi, real dz, real dx, real dt, int a, int nN, int nX, int nZ);
