#pragma once

#include <precision.hpp>
#include <variable.hpp>

inline int mod(int a, int b);
mode adamsBashforth(mode dfdt_current, mode dfdt_prev, real frac, real dt);
real checkCFL(Variable& psi, real dz, real dx, real dt, int a, int nN, int nX, int nZ);
