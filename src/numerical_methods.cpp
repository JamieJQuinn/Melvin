#include <cmath>
#include <iostream>
#include <numerical_methods.hpp>
#include <variable.hpp>

using std::cout;
using std::endl;
using std::isnan;

mode adamsBashforth(mode dfdt_current, mode dfdt_prev, real frac, real dt) {
  // Calcs X in equation T_{n+1} = T_{n} + X
  return ((2.0+frac)*dfdt_current - frac*dfdt_prev)*dt/2.0;
}

real checkCFL(Variable &psi, real dz, real dx, real dt, int a, int nN, int nX, int nZ) {
  real vxMax = 0.0f;
  real vzMax = 0.0f;
  real f=1.0f;

  psi.toPhysical();

  for(int k=0; k<nZ; ++k) {
    for(int j=0; j<nX; ++j) {
      real vx = psi.dfdzSpatial(j,k);
      real vz = psi.dfdx(j,k);
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

  if(dt > dz/vzMax or dt > dx/vxMax){
    cout << "CFL Condition Breached" << endl;
    exit(-1);
  }

  real threshold = 0.8;

  while(dt > threshold*dz/vzMax or dt > threshold*dx/vxMax) {
    dt*=threshold;
    f*=threshold;
  } 

  //if(vzMax < 0.2*dz/dt and vxMax < 0.2*(float(a)/nN)/dt) {
    //dt*=1.1;
    //f*=1.1;
  //} 

  if(f!=1.0f) {
    cout << "New time step is " << dt << endl;
  }
  return f;
}
