#pragma once

double dfdz(double *f, int k, double dz);
double dfdz2(double *f, int k, double dz);
double adamsBashforth(double dfdt_current, double dfdt_prev, double frac, double dt);
double checkCFL(double* psi, double dz, double dx, double dt, int a, int nN, int nX, int nZ);
