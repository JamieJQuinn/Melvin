#include <derivatives.hpp>
#include <cmath>

void computeLinearTemperatureDerivative(Variable &dTmpdt, const Variable &tmp, const Constants &c) {
  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      dTmpdt(n,k) = tmp.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*tmp(n,k);
    }
  }
}

void computeLinearVorticityDerivative(Variable &dOmgdt, const Variable &omg, const Variable &tmp, const Constants &c) {
  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      dOmgdt(n,k) =
        c.Pr*(omg.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*omg(n,k))
        + c.Pr*c.Ra*(n*M_PI/c.aspectRatio)*tmp(n,k);
    }
  }
}

void computeLinearXiDerivative(Variable &dXidt, const Variable &xi, Variable &dOmgdt, const Constants &c) {
  for(int n=0; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      dXidt(n,k) = c.tau*(xi.dfdz2(n,k) - pow(n*M_PI/c.aspectRatio, 2)*xi(n,k));
      dOmgdt(n,k) += -c.RaXi*c.tau*c.Pr*(n*M_PI/c.aspectRatio)*xi(n,k);
    }
  }
}

void computeLinearDerivatives(
    Variable &dTmpdt, const Variable &tmp,
    Variable &dOmgdt, const Variable &omg,
    Variable &dXidt, const Variable &xi,
    const Constants &c) {
  // Computes the (linear) derivatives of Tmp and omg
  computeLinearTemperatureDerivative(dTmpdt, tmp, c);
  computeLinearVorticityDerivative(dOmgdt, omg, tmp, c);
  if(c.isDoubleDiffusion) {
    computeLinearXiDerivative(dXidt, xi, dOmgdt, c);
  }
}

void addAdvectionApproximation(
    Variable &dTmpdt, const Variable &tmp,
    Variable &dOmgdt, const Variable &omg,
    Variable &dXidt, const Variable &xi,
    const Variable &psi,
    const Constants &c) {
  // Only applies to the linear simulation
  for(int k=1; k<c.nZ-1; ++k) {
    dOmgdt(0,k) = 0.0;
    dTmpdt(0,k) = 0.0;
  }
  for(int n=1; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      dTmpdt(n,k) += -1*tmp.dfdz(0,k)*n*M_PI/c.aspectRatio * psi(n,k);
    }
  }
  if(c.isDoubleDiffusion) {
    for(int k=1; k<c.nZ-1; ++k) {
      dXidt(0,k) = 0.0;
    }
    for(int n=1; n<c.nN; ++n) {
      for(int k=1; k<c.nZ-1; ++k) {
        dXidt(n,k) += -1*xi.dfdz(0,k)*n*M_PI/c.aspectRatio * psi(n,k);
      }
    }
  }
}

void computeNonlinearTemperatureDerivative(
    Variable &dTmpdt, const Variable &tmp,
    const Variable &psi,
    const Constants &c) {
  for(int n=1; n<c.nN; ++n) {
    for(int k=1; k<c.nZ-1; ++k) {
      // Contribution TO tmp[n=0]
      dTmpdt(0,k) +=
        -M_PI/(2*c.aspectRatio)*n*(
          psi.dfdz(n,k)*tmp(n,k) +
          tmp.dfdz(n,k)*psi(n,k)
          );
    }
  }
  #pragma omp parallel for schedule(dynamic)
  for(int n=1; n<c.nN; ++n) {
    // Contribution FROM tmp[n=0]
    for(int k=1; k<c.nZ-1; ++k) {
      dTmpdt(n,k) +=
        -n*M_PI/c.aspectRatio*psi(n,k)*tmp.dfdz(0,k);
    }
    // Contribution FROM tmp[n>0] and omg[n>0]
    int o;
    for(int m=1; m<n; ++m){
      // Case n = n' + n''
      o = n-m;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dTmpdt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          -m*psi.dfdz(o,k)*tmp(m,k)
          +o*tmp.dfdz(m,k)*psi(o,k)
          );
      }
    }
    for(int m=n+1; m<c.nN; ++m){
      // Case n = n' - n''
      o = m-n;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dTmpdt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          +m*psi.dfdz(o,k)*tmp(m,k)
          +o*tmp.dfdz(m,k)*psi(o,k)
          );
      }
    }
    for(int m=1; m+n<c.nN; ++m){
      // Case n= n'' - n'
      o = n+m;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dTmpdt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          +m*psi.dfdz(o,k)*tmp(m,k)
          +o*tmp.dfdz(m,k)*psi(o,k)
          );
      }
    }
  }
}

void computeNonlinearXiDerivative(
    Variable &dXidt, const Variable &xi,
    const Variable &psi,
    const Constants &c) {
  // Turns out it's the same for temperature and xi
  computeNonlinearTemperatureDerivative(dXidt, xi, psi, c);
}

void computeNonlinearVorticityDerivative(
    Variable &dOmgdt, const Variable &omg,
    const Variable &psi,
    const Constants &c) {
  #pragma omp parallel for schedule(dynamic)
  for(int n=1; n<c.nN; ++n) {
    int o;
    for(int m=1; m<n; ++m){
      // Case n = n' + n''
      o = n-m;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dOmgdt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          -m*psi.dfdz(o,k)*omg(m,k)
          +o*omg.dfdz(m,k)*psi(o,k)
          );
      }
    }
    for(int m=n+1; m<c.nN; ++m){
      // Case n = n' - n''
      o = m-n;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dOmgdt(n,k) +=
          -M_PI/(2*c.aspectRatio)*(
          +m*psi.dfdz(o,k)*omg(m,k)
          +o*omg.dfdz(m,k)*psi(o,k)
          );
      }
    }
    for(int m=1; m+n<c.nN; ++m){
      // Case n= n'' - n'
      o = n+m;
      assert(o>0 and o<c.nN);
      assert(m>0 and m<c.nN);
      for(int k=1; k<c.nZ-1; ++k) {
        dOmgdt(n,k) +=
          +M_PI/(2*c.aspectRatio)*(
          +m*psi.dfdz(o,k)*omg(m,k)
          +o*omg.dfdz(m,k)*psi(o,k)
          );
      }
    }
  }
}

void computeNonlinearDerivatives(
    Variable &dTmpdt, const Variable &tmp,
    Variable &dOmgdt, const Variable &omg,
    Variable &dXidt, const Variable &xi,
    const Variable &psi,
    const Constants &c) {
  computeNonlinearTemperatureDerivative(dTmpdt, tmp, psi, c);
  computeNonlinearVorticityDerivative(dOmgdt, omg, psi, c);
  if(c.isDoubleDiffusion) {
    computeNonlinearXiDerivative(dXidt, xi, psi, c);
  }
}
