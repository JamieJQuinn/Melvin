#include <iostream>
#include <cmath>
#include <cfloat>

const double EPSILON = FLT_EPSILON;

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
		// Constructor
		Sim(int Nz, int Nn, double dt,
		double Ra, double Pr, int a,
		double timeBetweenSaves, bool modifydt,
	       	int current, double t, double totalTime);

		// Helper functions
		double adamsBashforth(double *dfdt, int c, int Nz, int Nn,
			       	int k, int n, double dt, double frac);
		double dfdz2(double *f, int k, double oodz2);
		double triDiagonalSolver(const int Nz,
			       const double *rhs, double *sol, const double *sub,
			       const double * wk1, const double *wk2);
		double formTridiArrays ( const int Nz,
			const double *sub, const double *dia, const double *sup,
			double * wk1, double *wk2);

		void updateTmpAndOmg(double f);
		void computeLinearDerivatives();

		void runLinear();
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

	for(int i=0; i<Nz*Nn; ++i) {
		psi[i] = 0;
		omg[i] = 0;
		tmp[i] = 0;
	}
	for(int i=0; i<Nz*Nn*2; ++i) {
		dTmpdt[i] = 0;
		dOmgdt[i] = 0;
	}

	// Derived Constants
	Nx = Nz*a;
	dz = double(1)/(Nz-1);
	dx = double(a)/(Nx-1);
	oodz2 = pow(1.0/dz, 2);

}

double Sim::triDiagonalSolver(const int	Nz,
			       const double *rhs, double *sol, const double *sub,
			       const double * wk1, const double *wk2) {
	// Solves the tridiagonal system represented by sub, dia and sup.
	// If sub, dia and sup do not change, they can be rolled into wk1 and wk2
	// using formTridiArrays() and simply saved

	// Forward Subsitution
	sol[0] = rhs[0]*wk1[0];
	for (int i=1; i<Nz; ++i) {
		sol[i] = (rhs[i] - sub[i-1]*sol[i-1])*wk1[i];
	}
	// Backward Substitution
	for (int i=Nz-2; i>=0; --i) {
		sol[i] -= wk2[i]*sol[i+1];
	}
}

double Sim::formTridiArrays ( const int Nz,
	       	const double *sub, const double *dia, const double *sup,
		double * wk1, double *wk2) {
	wk1[0] = 1.0/dia[0];
	wk2[0] = sup[0]*wk1[0];
	for (int i=1; i<Nz-1; ++i) {
		wk1[i] = 1.0/(dia[i] - sub[i-1] * wk2[i-1]);
		wk2[i] = sup[i]*wk1[i];
	}
	wk1[Nz-1] = 1.0/(dia[Nz-1] - sub[Nz-2]*wk2[Nz-2]);
	// no value for wk2[Nz-1], it never gets called anyway
}

bool test_Sim_triDiagonalSolver() {
	Sim sim = Sim(101, 51, 1.0e-6, 1e6, 0.5, 3, 1.5e-3, true, 0, 0, 1);
	double rhs[] = {3.0, 5.0, 3.0};
	double sol[3];
	double sub[] = {2.0, 2.0};
	double dia[] = {1.0, 1.0, 1.0};
	double wk1[3];
	double wk2[3];
	sim.formTridiArrays(3, sub, dia, sub, wk1, wk2);
	sim.triDiagonalSolver(3, rhs, sol, sub, wk1, wk2);
	bool pass1 = true;
	for(int i=0; i<3; ++i) {
		if(sol[i] - 1.0 > EPSILON) {
			pass1 = false;
		}
	}
	rhs[0] = rhs[2] = -1.0;
	rhs[1] = 0.0;
	sub[0] = sub[1] = 1.0;
	dia[0] = dia[1] = dia[2] = -2.0;
	sim.formTridiArrays(3, sub, dia, sub, wk1, wk2);
	sim.triDiagonalSolver(3, rhs, sol, sub, wk1, wk2);
	bool pass2 = true;
	for(int i=0; i<3; ++i) {
		if(sol[i] - 1.0 > EPSILON) {
			pass2 = false;
		}
	}
	return pass1 and pass2;
}

double Sim::adamsBashforth(double *dfdt, int c, int Nz, int Nn, int k, int n, double dt, double frac) {
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

double Sim::dfdz2(double *f, int k, double oodz2) {
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
	}	double sol[3];


}

int main() {
	Sim simulation = Sim(101, 51, 1.0e-6, 1e6, 0.5, 3, 1.5e-3, true, 0, 0, 1);
	simulation.runLinear();

	// cout << test_Sim_triDiagonalSolver() << endl;
}

