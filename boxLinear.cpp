#include <iostream>
#include <cmath>
using namespace std;

class Sim {
	public:
		// Defined constants
		int Nz;
		int Nn; 
		double dt;
		double Ra; 
		double Pr;
	        int a;
		double timeBetweenSaves; 
		bool modifydt;
	       	int current; 
		double t; 
		double totalTime;
		// Derived constants
		double dz;
		double dx;
		int Nx;
		double oodz2;
		// Variable arrays
		double * psi; // Stream function (Psi)
		double * omg; // Vorticity (Omega)
		double * tmp; // Temperature

		double * dOmgdt; // d/dt of vorticty
		double * dTmpdt; // d/dt of temperature

		Sim(int Nz, int Nn, double dt,
		double Ra, double Pr, int a,
		double timeBetweenSaves, bool modifydt,
	       	int current, double t, double totalTime);

		void runLinear();
		void updateTmpAndOmg(double f);
		void computeLinearDerivatives();
};

Sim::Sim(int Nz, int Nn, double dt,
		double Ra, double Pr, int a,
		double timeBetweenSaves, bool modifydt,
	       	int current, double t, double totalTime) 
	: Nz {Nz}
	, Nn {Nn}
	, dt {dt}
	, Ra {Ra}
	, a {a}
	, timeBetweenSaves {timeBetweenSaves}
	, modifydt {modifydt}
	, current {current}
	, t {t}
	, totalTime {totalTime}
{
	psi = new double [Nn*Nz];
	omg = new double [Nn*Nz];
	tmp = new double [Nn*Nz];

	dTmpdt = new double [2*Nn*Nz];
	dOmgdt = new double [2*Nn*Nz];

	Nx = Nz*a;
	dz = double(1)/(Nz-1);
	dx = double(a)/(Nx-1);
	oodz2 = pow(1.0/dz, 2);


	for(int i=0; i<Nz*Nn; ++i) {
		psi[i] = 0;
		omg[i] = 0;
		tmp[i] = 0;
	}
	for(int i=0; i<Nz*Nn*2; ++i) {
		dTmpdt[i] = 0;
		dOmgdt[i] = 0;
	}

}
	

double adamsBashforth(double *dfdt, int c, int Nz, int Nn, int k, int n, double dt, double frac) {
	return ((1+frac/2)*dfdt[c*Nz*Nn+n*Nz+k] - frac/2*dfdt[((c+1)%2)*Nz*Nn+n*Nz+k])*dt;
}


void Sim::updateTmpAndOmg(double f = 1.0) {
	for(int n=0; n<Nn; ++n) {
		for(int k=0; k<Nz; ++k) {
			tmp[n*Nz+k] = adamsBashforth(dTmpdt, current, Nz, Nn, k, n, dt, f);
			omg[n*Nz+k] = adamsBashforth(dOmgdt, current, Nz, Nn, k, n, dt, f);
		}
	}
}

double dfdz2(double *f, int k, double oodz2) {
	return (f[k+1] - 2*f[k] + f[k+1])*oodz2;
}

void Sim::computeLinearDerivatives() {
	for(int n=0; n<Nn; ++n) {
		for(int k=1; k<Nz-1; ++k) {
			int di = current*Nz*Nn+n*Nz+k;
			int i = k+n*Nz;
			dTmpdt[di] = dfdz2(tmp, i, oodz2)
				- pow(n*M_PI/a, 2)*tmp[i]
				+ n*M_PI/a * psi[i];
			dOmgdt[di] = Pr*(dfdz2(omg, i, oodz2)
				- pow(n*M_PI/a, 2)*omg[i])
				+ Ra*Pr*n*M_PI/a*tmp[i];
		}
	}
}

void Sim::runLinear() {
	// Initial Conditions
	// Let psi = omg = dtmpdt = domgdt = 0
	// Let tmp[n=1] = sin(PI*z)
	// and tmp[n=0] = (1-z)/N
	for(int k=0; k<Nz; ++k) {
		tmp[0] = 1-k*dz;
		tmp[1] = sin(M_PI*k*dz);
	}
}

int main() {
	Sim simulation = Sim(101, 51, 1.0e-6, 1e6, 0.5, 3, 1.5e-3, true, 0, 0, 1);
	simulation.runLinear();
}

