#include <cmath>
#include "numerical_methods.hpp"

double adamsBashforth(double dfdt_current, double dfdt_prev, double frac, double dt) {
  // Calcs X in equation T_{n+1} = T_{n} + X
  return ((2.0+frac)*dfdt_current - frac*dfdt_prev)*dt/2.0;
}

double dfdz(double *f, int k, double dz) {
  return (f[k+1]-f[k-1])/(2*dz);
}

double dfdz2(double *f, int k, double dz) {
  return (f[k+1] - 2*f[k] + f[k-1])/pow(dz,2);
}

