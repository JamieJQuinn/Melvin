#pragma once

#include <precision.hpp>
#include <variable.hpp>

inline int mod(int a, int b);
real adamsBashforth(real dfdt_current, real dfdt_prev, real frac, real dt);
real checkCFL(const Variable& psi, real dz, real dx, real dt, int a, int nN, int nX, int nZ);
