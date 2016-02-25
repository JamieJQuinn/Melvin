#include <iostream>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <string>
#include <cassert>

const double EPSILON = DBL_EPSILON;
#define DEBUG

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

		double * wk1; // Used in Thomas Algorithm
		double * wk2;
		double * sub;

		// Constructor
		Sim(int Nz, int Nn, double dt,
		double Ra, double Pr, int a,
		double timeBetweenSaves, bool modifydt,
	       	int current, double t, double totalTime);

		// Destructor
		~Sim();

		// Helper functions
		double adamsBashforth(double *dfdt,int k, int n, double frac);
		double dfdz2(double *f, int k);
		double triDiagonalSolver(const int Nz,
			       const double *rhs, double *sol, const double *sub,
			       const double * wk1, const double *wk2);
		double formTridiArrays ( const int Nz,
			const double *sub, const double *dia, const double *sup,
			double * wk1, double *wk2);
		void printMaxOf(double *a, std::string name);

		// Simulation functions
		void updateTmpAndOmg(double f);
		void computeLinearDerivatives();
		void solveForPsi();

		// Runs the linear simulation
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
	, Pr {Pr}
	, a {a}
	, timeBetweenSaves {timeBetweenSaves}
	, modifydt {modifydt}
	, current {current}
	, t {t}
	, totalTime {totalTime}
{
	// Derived Constants
	Nx = Nz*a;
	dz = double(1)/(Nz-1);
	dx = double(a)/(Nx-1);
	oodz2 = pow(1.0/dz, 2);

	// Initialise Arrays
	psi = new double [Nn*Nz];
	omg = new double [Nn*Nz];
	tmp = new double [Nn*Nz];

	dTmpdt = new double [2*Nn*Nz];
	dOmgdt = new double [2*Nn*Nz];

	wk1 = new double [Nn*Nz];
	wk2 = new double [Nn*Nz];
	sub = new double [Nz];

	for(int i=0; i<Nz*Nn; ++i) {
		psi[i] = 0.0;
		omg[i] = 0.0;
		tmp[i] = 0.0;
	}
	for(int i=0; i<Nz*Nn*2; ++i) {
		dTmpdt[i] = 0.0;
		dOmgdt[i] = 0.0;
	}

	// Precalculate tridiagonal stuff
	double * dia = new double [Nz];
	double * sup = new double [Nz];
	for(int k=0; k<Nz; ++k) {
		sub[k] = sup[k] = -oodz2;
	}
	for(int n=0; n<Nn; ++n) {
		for(int k=0; k<Nz; ++k){
			dia[k] = pow(M_PI/a*n, 2) + 2*oodz2;
		}
		dia[0] = dia[Nz-1] = 1.0;
		sub[Nz-2] = sup[0] = 0.0;
		formTridiArrays( Nz,
		sub, dia, sup,
		wk1+Nz*n, wk2+n*Nz);
	}
	for(int i=0; i<Nz*Nn; ++i) {
		assert(!std::isnan(wk1[i]));
		assert(!std::isnan(wk2[i]));
	}
	delete [] dia;
	delete [] sup;
}

Sim::~Sim() {
	// Destructor
	delete[] psi  ;
	delete[] omg  ;
	delete[] tmp  ;

	delete[] dTmpdt  ;
	delete[] dOmgdt  ;

	delete[] wk1  ;
	delete[] wk2  ;
	delete[] sub;
}

double Sim::triDiagonalSolver(const int	Nz,
			       const double *rhs, double *sol, const double *sub,
			       const double * wk1, const double *wk2) {
	// Solves the tridiagonal system represented by sub, dia and sup.
	// If sub, dia and sup do not change, they can be rolled into wk1 and wk2
	// using formTridiArrays() and simply saved

	// Forward Subsitution
	assert(!isnan(wk1[0]));
	assert(!isnan(rhs[0]));
	sol[0] = rhs[0]*wk1[0];
	for (int i=1; i<Nz; ++i) {
		sol[i] = (rhs[i] - sub[i-1]*sol[i-1])*wk1[i];
		assert(!isnan(rhs[i]));
		assert(!isnan(sub[i-1]));
		assert(!isnan(sol[i-1]));
		assert(!isnan(wk1[i]));
		assert(!isnan(sol[i]));
	}
	// Backward Substitution
	for (int i=Nz-2; i>=0; --i) {
		sol[i] -= wk2[i]*sol[i+1];
	}
}

double Sim::formTridiArrays ( const int Nz,
	       	const double *sub, const double *dia, const double *sup,
		double * wk1, double *wk2) {
	assert(dia[0] != 0.0);
	wk1[0] = 1.0/dia[0];
	wk2[0] = sup[0]*wk1[0];
	for (int i=1; i<Nz-1; ++i) {
		assert((dia[i] - sub[i-1] * wk2[i-1]) != 0.0);
		wk1[i] = 1.0/(dia[i] - sub[i-1] * wk2[i-1]);
		wk2[i] = sup[i]*wk1[i];
	}
	assert((dia[Nz-1] - sub[Nz-2]*wk2[Nz-2]) != 0.0);
	wk1[Nz-1] = 1.0/(dia[Nz-1] - sub[Nz-2]*wk2[Nz-2]);
}

bool test_Sim_triDiagonalSolver() {
	Sim sim = Sim(101, 51, 1.0e-6, 1e6, 0.5, 3, 1.5e-3, true, 0, 0, 1);
	// Test 1
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
	// Test 2
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

double Sim::adamsBashforth(double *dfdt, int k, int n, double frac) {
	// Calcs X in equation T_{n+1} = T_{n} + X
	return ((2.0+frac)*dfdt[current*Nz*Nn+n*Nz+k] - frac*dfdt[((current+1)%2)*Nz*Nn+n*Nz+k])*dt/2.0;
}

void Sim::updateTmpAndOmg(double f = 1.0) {
	// Update variables using Adams-Bashforth Scheme
	// f is the proportional change between the new dt and old dt
	// ( if dt changed )
	for(int n=0; n<Nn; ++n) {
		for(int k=0; k<Nz; ++k) {
			tmp[n*Nz+k] += adamsBashforth(dTmpdt, k, n, f);
			omg[n*Nz+k] += adamsBashforth(dOmgdt, k, n, f);
		}
		// check BCs
		if(n>0) {
			assert(tmp[n*Nz] < EPSILON);
		} else {
			assert(tmp[n*Nz] - 1.0 < EPSILON);
		}
		assert(tmp[n*Nz+Nz-1] < EPSILON);
		assert(omg[n*Nz] < EPSILON);
		assert(omg[n*Nz+Nz-1] < EPSILON);
	}	
}

double Sim::dfdz2(double *f, int k) {
	// Calculates d/dz of f
	return (f[k+1] - 2*f[k] + f[k-1])*oodz2;
}

bool test_Sim_dfdz2() {
	int NTests = 1;
	bool pass [NTests];
	for(int i = 0; i<NTests; ++i){
		pass[i] = true;
	}
	Sim sim = Sim(101, 51, 1.0e-6, 1e6, 0.5, 3, 1.5e-3, true, 0, 0, 1);
	double T [] = {3.0, 2.0, 3.0};
	
	if (sim.dfdz2(T, 1)/sim.oodz2 - 2.0 > EPSILON) {
		pass[0] = false;
	}

	bool passTotal = true;
	for(int i=0; i<NTests; ++i) {
		passTotal &= pass[i];
	}
	return passTotal;
}

void Sim::computeLinearDerivatives() {
	// Computes the (linear) derivatives of Tmp and omg for n>0
	for(int n=1; n<Nn; ++n) {
		for(int k=1; k<Nz-1; ++k) {
			// Setup indices
			int di = current*Nz*Nn+n*Nz+k;
			int i = k+n*Nz;

			dTmpdt[di] = 
				dfdz2(tmp, i) - pow(n*M_PI/a, 2)*tmp[i] 
				+ n*M_PI/a * psi[i];
			dOmgdt[di] = 
				Pr*(
				dfdz2(omg, i) 
				- pow(n*M_PI/a, 2)*omg[i] 
				+ Ra*n*M_PI/a*tmp[i]
				);
		}
	}
}

void Sim::solveForPsi(){
	// Solves for Psi using Thomas algorithm
	for(int n=0; n<Nn; ++n) {
		triDiagonalSolver(Nz, omg+Nz*n, psi+Nz*n, sub, wk1+Nz*n, wk2+Nz*n);
		// Check Boundary Conditions
		assert(psi[Nz*n+0] == 0.0);
		assert(psi[Nz*n+Nz-1] == 0.0);
	}
	for(int k=0; k<Nz; ++k) {
		// Check BCs
		assert(psi[Nz*0+k] < EPSILON);
	}

}

void Sim::printMaxOf(double *a, std::string name) {
	int nStart = 1; // n level to start from
	// Find max
	double max = a[nStart*Nz];
	int maxLoc[2];
	for(int n=nStart; n<Nn; ++n) {
		for(int k=0; k<Nz; ++k) {
			if(a[n*Nz+k]>max) {
				max = a[n*Nz+k];
				maxLoc[0] = k;
				maxLoc[1] = n;
			}
		}
	}
	// print max
	printf("%e @ (%d, %d)", max, maxLoc[0], maxLoc[1]);
}


void Sim::runLinear() {
	// Initial Conditions
	// Let psi = omg = dtmpdt = domgdt = 0
	// Let tmp[n>0] = sin(PI*z)
	// and tmp[n=0] = (1-z)/N
	for(int k=0; k<Nz; ++k) {
		tmp[Nz*0+k] = 1-k*dz;
		for(int n=1; n<Nn; ++n) {
			tmp[Nz*n+k] = sin(M_PI*k*dz);
		}
	}	
	// Check BCs
	for(int n=0; n<Nn; ++n){
		if(n>0) {
			assert(tmp[n*Nz] < EPSILON);
		} else {
			assert(tmp[n*Nz] - 1.0 < EPSILON);
		}
		assert(tmp[n*Nz+Nz-1] < EPSILON);
		assert(omg[n*Nz] < EPSILON);
		assert(omg[n*Nz+Nz-1] < EPSILON);
	}

	// Stuff for critical rayleigh check
	double tmpPrev[Nn];
	double omgPrev[Nn];
	double psiPrev[Nn];
	for(int n=0; n<Nn; ++n){
		tmpPrev[n] = tmp[32+n*Nz];
		psiPrev[n] = psi[32+n*Nz];
		omgPrev[n] = omg[32+n*Nz];
	}
	current = 0;
	int steps = 0;
	while (t<totalTime) {
		if(steps%500 == 0) {
			for(int n=1; n<Nn; ++n){
				printf("%d: %e, %e, %e\n",n,
						std::log(std::abs(tmp[32+n*Nz])) - std::log(std::abs(tmpPrev[n])),
						std::log(std::abs(omg[32+n*Nz])) - std::log(std::abs(omgPrev[n])),
						std::log(std::abs(psi[32+n*Nz])) - std::log(std::abs(psiPrev[n])));
				tmpPrev[n] = tmp[32+n*Nz];
				psiPrev[n] = psi[32+n*Nz];
				omgPrev[n] = omg[32+n*Nz];
				/*
				printf("%d: %e, %e, %e\n", n, 
						tmp[Nz*n+32],
						omg[Nz*n+32],
						psi[Nz*n+32]);
						*/

			}
			/*
			printMaxOf(tmp, "tmp");
			printMaxOf(omg, "omg");
			printMaxOf(psi, "psi");
			//printMaxOf(dOmgdt+current*Nn*Nz, "dOmgdt");
			std::printf(" \n");
			*/
		}
		steps++;
		computeLinearDerivatives();
		updateTmpAndOmg();
		solveForPsi();
		t+=dt;
		++current%=2;
	}	
}

int main() {
// Sim::Sim(int Nz, int Nn, double dt, double Ra, double Pr, int a,double timeBetweenSaves, bool modifydt, int current, double t, double totalTime
	Sim simulation = Sim(101, 10, 1.0e-5, 660.5, 1, 3, 1.5e-3, false, 0, 0, 1e1);
	simulation.runLinear();

	// test_Sim_triDiagonalSolver();
	// test_Sim_dfdz2();
}

