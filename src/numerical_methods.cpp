#include <cmath>
#include <iostream>
#include "numerical_methods.hpp"

using namespace std;

real adamsBashforth(real dfdt_current, real dfdt_prev, real frac, real dt) {
  // Calcs X in equation T_{n+1} = T_{n} + X
  return ((2.0+frac)*dfdt_current - frac*dfdt_prev)*dt/2.0;
}

real dfdz(real *f, int k, real dz) {
  return (f[k+1]-f[k-1])/(2*dz);
}

real dfdz2(real *f, int k, real dz) {
  return (f[k+1] - 2*f[k] + f[k-1])/pow(dz,2);
}

real checkCFL(real* psi, real dz, real dx, real dt, int a, int nN, int nX, int nZ) {
  real vxMax = 0.0f;
  real vzMax = 0.0f;
  real f=1.0f;
  for(int j=1; j<nX-1; ++j) {
    for(int k=1; k<nZ-1; ++k) {
      real vx = 0.0;
      real vz = 0.0;
      for(int n=0; n<nN; ++n) {
        vx += dfdz(psi, n*nZ+k, dz)*sin(n*M_PI*j*dx/a);
        vz += n*M_PI/a*psi[n*nZ+k]*cos(n*M_PI*j*dx/a);
      }
      if(isnan(vx) or isnan(vz)){
        cout << "CFL Condition Breached" << endl;
        exit(-1);
      }
      if( std::abs(vx) > vxMax ) {
        vxMax = std::abs(vx);
      }
      if( std::abs(vz) > vzMax ) {
        vzMax = std::abs(vz);
      }
    }
  }
  if(vzMax > dz/dt or vxMax > (float(a)/nN)/dt){
    cout << "CFL Condition Breached" << endl;
    exit(-1);
  } 
  while(vzMax > 0.9*dz/dt or vxMax > 0.9*(float(a)/nN)/dt) {
    dt*=0.9;
    f*=0.9;
  } 
  if(f!=1.0f) {
    cout << "New time step is " << dt << endl;
  }
  return f;
}
