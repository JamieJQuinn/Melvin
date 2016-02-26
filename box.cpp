#include <iostream>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <string>
#include <cassert>

const double EPSILON = DBL_EPSILON;
const int THREAD_COUNT = 8;
#define NDEBUG

class Sim {
	public:
		// Defined constants
		int nZ;
		int nN; 
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
		Sim(int nZ, int nN, double dt,
		double Ra, double Pr, int a,
		double timeBetweenSaves, bool modifydt,
	       	int current, double t, double totalTime);

		// Destructor
		~Sim();

		// Helper functions
		double adamsBashforth(double *dfdt,int k, int n, double frac);
		double inline dfdz(double *f, int k);
		double inline dfdz2(double *f, int k);
		double triDiagonalSolver(const int nZ,
			       const double *rhs, double *sol, const double *sub,
			       const double * wk1, const double *wk2);
		double formTridiArrays ( const int nZ,
			const double *sub, const double *dia, const double *sup,
			double * wk1, double *wk2);
		void printMaxOf(double *a, std::string name);
		void printBenchmarkData();

		// Simulation functions
		void updateTmpAndOmg(double f);
		void computeLinearDerivatives(int linearSim = 1);
		void computeNonLinearDerivatives();
		void computeNonLinearDerivativesFor(int n);
		void solveForPsi();

		// Runs the linear simulation
		void runLinear();
		void runNonLinear();
};

Sim::Sim(int nZ, int nN, double dt,
		double Ra, double Pr, int a,
		double timeBetweenSaves, bool modifydt,
	       	int current, double t, double totalTime) 
	: nZ {nZ}
	, nN {nN}
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
	Nx = nZ*a;
	dz = double(1)/(nZ-1);
	dx = double(a)/(Nx-1);
	oodz2 = pow(1.0/dz, 2);

	// Initialise Arrays
	psi = new double [nN*nZ];
	omg = new double [nN*nZ];
	tmp = new double [nN*nZ];

	dTmpdt = new double [2*nN*nZ];
	dOmgdt = new double [2*nN*nZ];

	wk1 = new double [nN*nZ];
	wk2 = new double [nN*nZ];
	sub = new double [nZ];

	for(int i=0; i<nZ*nN; ++i) {
		psi[i] = 0.0;
		omg[i] = 0.0;
		tmp[i] = 0.0;
	}
	for(int i=0; i<nZ*nN*2; ++i) {
		dTmpdt[i] = 0.0;
		dOmgdt[i] = 0.0;
	}

	// Precalculate tridiagonal stuff
	double * dia = new double [nZ];
	double * sup = new double [nZ];
	for(int k=0; k<nZ; ++k) {
		sub[k] = sup[k] = -oodz2;
	}
	for(int n=0; n<nN; ++n) {
		for(int k=0; k<nZ; ++k){
			dia[k] = pow(M_PI/a*n, 2) + 2*oodz2;
		}
		dia[0] = dia[nZ-1] = 1.0;
		sub[nZ-2] = sup[0] = 0.0;
		formTridiArrays( nZ,
		sub, dia, sup,
		wk1+nZ*n, wk2+n*nZ);
	}
	for(int i=0; i<nZ*nN; ++i) {
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

double Sim::triDiagonalSolver(const int	nZ,
			       const double *rhs, double *sol, const double *sub,
			       const double * wk1, const double *wk2) {
	// Solves the tridiagonal system represented by sub, dia and sup.
	// If sub, dia and sup do not change, they can be rolled into wk1 and wk2
	// using formTridiArrays() and simply saved

	// Forward Subsitution
	assert(!isnan(wk1[0]));
	assert(!isnan(rhs[0]));
	sol[0] = rhs[0]*wk1[0];
	for (int i=1; i<nZ; ++i) {
		sol[i] = (rhs[i] - sub[i-1]*sol[i-1])*wk1[i];
		assert(!isnan(rhs[i]));
		assert(!isnan(sub[i-1]));
		assert(!isnan(sol[i-1]));
		assert(!isnan(wk1[i]));
		assert(!isnan(sol[i]));
	}
	// Backward Substitution
	for (int i=nZ-2; i>=0; --i) {
		sol[i] -= wk2[i]*sol[i+1];
	}
}

double Sim::formTridiArrays ( const int nZ,
	       	const double *sub, const double *dia, const double *sup,
		double * wk1, double *wk2) {
	assert(dia[0] != 0.0);
	wk1[0] = 1.0/dia[0];
	wk2[0] = sup[0]*wk1[0];
	for (int i=1; i<nZ-1; ++i) {
		assert((dia[i] - sub[i-1] * wk2[i-1]) != 0.0);
		wk1[i] = 1.0/(dia[i] - sub[i-1] * wk2[i-1]);
		wk2[i] = sup[i]*wk1[i];
	}
	assert((dia[nZ-1] - sub[nZ-2]*wk2[nZ-2]) != 0.0);
	wk1[nZ-1] = 1.0/(dia[nZ-1] - sub[nZ-2]*wk2[nZ-2]);
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
	return ((2.0+frac)*dfdt[current*nZ*nN+n*nZ+k] - frac*dfdt[((current+1)%2)*nZ*nN+n*nZ+k])*dt/2.0;
}

void Sim::updateTmpAndOmg(double f = 1.0) {
	// Update variables using Adams-Bashforth Scheme
	// f is the proportional change between the new dt and old dt
	// ( if dt changed )
	for(int n=0; n<nN; ++n) {
		for(int k=0; k<nZ; ++k) {
			tmp[n*nZ+k] += adamsBashforth(dTmpdt, k, n, f);
			omg[n*nZ+k] += adamsBashforth(dOmgdt, k, n, f);

			assert(!isnan(tmp[n*nZ+k]));
			assert(!isnan(omg[n*nZ+k]));
		}
		// check BCs
		if(n>0) {
			assert(tmp[n*nZ] < EPSILON);
		} else {
			assert(tmp[n*nZ] - 1.0 < EPSILON);
		}
		assert(tmp[n*nZ+nZ-1] < EPSILON);
		assert(omg[n*nZ] < EPSILON);
		assert(omg[n*nZ+nZ-1] < EPSILON);
	}	
}

double Sim::dfdz(double *f, int k) {
	assert(!isnan(f[k+1]));
	assert(!isnan(f[k-1]));
	assert(!isnan(1.0/(2*dz)));
	assert(!isnan((f[k+1]-f[k-1])/(2*dz)));
	return (f[k+1]-f[k-1])/(2*dz);
}

bool test_Sim_dfdz() {
	int NTests = 2;
	bool pass [NTests];
	for(int i = 0; i<NTests; ++i){
		pass[i] = true;
	}

	Sim sim = Sim(101, 51, 1.0e-6, 1e6, 0.5, 3, 1.5e-3, true, 0, 0, 1);
	double T [] = {3.0, 2.0, 3.0};
	
	if (sim.dfdz(T, 1)*(2*sim.dz) - 0.0 > EPSILON) {
		pass[0] = false;
	}

	T[2] = 1.0;
	if (sim.dfdz(T, 1)*(2*sim.dz) - (-2) > EPSILON) {
		pass[1] = false;
	}

	bool passTotal = true;
	for(int i=0; i<NTests; ++i) {
		passTotal &= pass[i];
	}
	return passTotal;
}

double Sim::dfdz2(double *f, int k) {
	assert(!isnan((f[k+1] - 2*f[k] + f[k-1])*oodz2));
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

void Sim::computeLinearDerivatives(int linearSim) {
	// Computes the (linear) derivatives of Tmp and omg
	// If linear sim is 0, we start n from 0 and the advection approximation
	// in dTmpdt vanishes
	for(int n=linearSim; n<nN; ++n) {
		for(int k=1; k<nZ-1; ++k) {
			// Setup indices
			int di = current*nZ*nN+n*nZ+k;;

			int i = k+n*nZ;

			dTmpdt[di] = dfdz2(tmp, i) - pow(n*M_PI/a, 2)*tmp[i];
			if (linearSim == 1) {
				dTmpdt[di] += n*M_PI/a * psi[i];
			}
			assert(!isnan(dfdz2(tmp, i) - pow(n*M_PI/a, 2)*tmp[i] 
				+ n*M_PI/a * psi[i]*linearSim));
			dOmgdt[di] = 
				Pr*(dfdz2(omg, i) - pow(n*M_PI/a, 2)*omg[i] 
				+ Ra*n*M_PI/a*tmp[i]
				);
			assert(!isnan(Pr*(
				dfdz2(omg, i) 
				- pow(n*M_PI/a, 2)*omg[i] 
				+ Ra*n*M_PI/a*tmp[i]
				)));
			assert(dOmgdt[nZ*0+k] < EPSILON);
		}
	}
}

void Sim::computeNonLinearDerivativesFor(int n){
		for(int k=1; k<nZ-1; ++k) {
			int in = n*nZ + k;
			// Contribution TO tmp[n=0]
			dTmpdt[current*nZ*nN+0*nN+k] += 
				-M_PI/(2*a)*n*(
					dfdz(psi, in)*tmp[in] +
					dfdz(tmp, in)*psi[in]
					);
			// Contribution FROM tmp[n=0]
			dTmpdt[current*nZ*nN + in] += 
				-n*M_PI/a*psi[in]*dfdz(tmp, 0*nZ+k);
		}
		// Contribution FROM tmp[n>0] and omg[n>0]
		for(int m=1; m<nN; ++m){
			for(int k=1; k<nZ-1; ++k) {
				int im = nZ*m+k;
				// Case n = n' + n''
				int o = n-m; 
				int io = nZ*o + k;
				if(o > 0 && o < nN) {
					dTmpdt[current*nZ*nN+nZ*n+k] += 
						-M_PI/(2*a)*(
						-m*dfdz(psi, io)*tmp[im]
						+o*dfdz(tmp, im)*psi[io]
						);
					dOmgdt[current*nZ*nN+nZ*n+k] += 
						-M_PI/(2*a)*(
						-m*dfdz(psi, io)*omg[im]
						+o*dfdz(omg, im)*psi[io]
						);
				}
				// Case n = n' - n''
				o = m-n;
				io = nZ*o + k;
				if(o > 0 && o < nN) {
					dTmpdt[current*nZ*nN+nZ*n+k] += 
						-M_PI/(2*a)*(
						+m*dfdz(psi, io)*tmp[im]
						+o*dfdz(tmp, im)*psi[io]
						);
					dOmgdt[current*nZ*nN+nZ*n+k] += 
						-M_PI/(2*a)*(
						+m*dfdz(psi, io)*omg[im]
						+o*dfdz(omg, im)*psi[io]
						);
				}
				// Case n= n'' - n'
				o = n+m;
				io = nZ*o + k;
				if(o > 0 && o < nN) {
					dTmpdt[current*nZ*nN+nZ*n+k] += 
						-M_PI/(2*a)*(
						+m*dfdz(psi, io)*tmp[im]
						+o*dfdz(tmp, im)*psi[io]
						);
					dOmgdt[current*nZ*nN+nZ*n+k] += 
						+M_PI/(2*a)*(
						+m*dfdz(psi, io)*omg[im]
						+o*dfdz(omg, im)*psi[io]
						);
				}
			}
		}
}

void Sim::computeNonLinearDerivatives() { 
	#pragma omp parallel for
	for(int n=1; n<nN; ++n) {
		computeNonLinearDerivativesFor(n);
	}
}

void Sim::solveForPsi(){
	// Solve for Psi using Thomas algorithm
	for(int n=0; n<nN; ++n) {
		triDiagonalSolver(nZ, omg+nZ*n, psi+nZ*n, sub, wk1+nZ*n, wk2+nZ*n);
		// Check Boundary Conditions
		assert(psi[nZ*n+0] == 0.0);
		assert(psi[nZ*n+nZ-1] == 0.0);
	}
	// Check BCs
	for(int k=0; k<nZ; ++k) {
		assert(psi[nZ*0+k] < EPSILON);
	}

}

void Sim::printMaxOf(double *a, std::string name) {
	int nStart = 0; // n level to start from
	// Find max
	double max = a[nStart*nZ];
	int maxLoc[] = {0, nStart};
	for(int n=nStart; n<nN; ++n) {
		for(int k=0; k<nZ; ++k) {
			if(a[n*nZ+k]>max) {
				max = a[n*nZ+k];
				maxLoc[0] = k;
				maxLoc[1] = n;
			}
		}
	}
	// print max
	printf("%e @ (%d, %d)", max, maxLoc[0], maxLoc[1]);
}

void Sim::printBenchmarkData() {
	for(int n=0; n<21; ++n) {
		printf("%d | %e | %e | %e\n", n, tmp[n*nZ+33], omg[n*nZ+33], psi[n*nZ+33]);
	}
}

void Sim::runNonLinear() {
	// Initial Conditions
	// Let psi = omg = dtmpdt = domgdt = 0
	// Let tmp[n] = 0.01*sin(PI*z) for certain n
	// and tmp[n=0] = (1-z)/N
	int nInit [] = {1}; 
	for(int k=0; k<nZ; ++k) {
		tmp[nZ*0+k] = 1-k*dz;
		for(int n:nInit) {
			tmp[nZ*n+k] = 0.01f*sin(M_PI*k*dz);
		}
	}
	current = 0;
	int steps = 0;
	while (t<totalTime) {
		//printf("%e\n", t);
		if(steps%1000 == 0) {
			for(int n=1; n<nN; ++n){
				/*
				printf("%d: %e, %e, %e\n", n, 
						tmp[nZ*n+32],
						omg[nZ*n+32],
						psi[nZ*n+32]);
						*/

			}
			/*
			int points = 5;
			int kSpots [] = {17, 33, 50, 66, 84};
			int nSpots [] = {14};
			for(int n: nSpots) {
				for(int k:kSpots){
					printf("%e|", tmp[k+n*nZ]);
				}
				printf("\n");
			}
			*/
			/*
			n=9;
			printf("%e|%e|%e|%e|%e|%e\n", tmp[n*nZ+0],tmp[n*nZ+(nZ-1)/5],tmp[n*nZ+2*(nZ-1)/5],tmp[n*nZ+3*(nZ-1)/5],tmp[n*nZ+4*(nZ-1)/5],tmp[n*nZ+nZ-1]);
			*/
			/*
			printMaxOf(tmp, "tmp");
			printMaxOf(omg, "omg");
			printMaxOf(psi, "psi");
			//printMaxOf(dOmgdt+current*nN*nZ, "dOmgdt");
			//printMaxOf(dTmpdt+current*nN*nZ, "dOmgdt");
			std::printf(" \n");
			*/
			printBenchmarkData();
			std::cout << std::endl;
		}
		steps++;
		computeLinearDerivatives(0);
		computeNonLinearDerivatives();
		updateTmpAndOmg();
		solveForPsi();
		t+=dt;
		++current%=2;
	}	
	printBenchmarkData();

}

void Sim::runLinear() {
	// Initial Conditions
	// Let psi = omg = dtmpdt = domgdt = 0
	// Let tmp[n>0] = sin(PI*z)
	// and tmp[n=0] = (1-z)/N
	for(int k=0; k<nZ; ++k) {
		tmp[nZ*0+k] = 1-k*dz;
		for(int n=1; n<nN; ++n) {
			tmp[nZ*n+k] = sin(M_PI*k*dz);
		}
	}	
	// Check BCs
	for(int n=0; n<nN; ++n){
		if(n>0) {
			assert(tmp[n*nZ] < EPSILON);
		} else {
			assert(tmp[n*nZ] - 1.0 < EPSILON);
		}
		assert(tmp[n*nZ+nZ-1] < EPSILON);
		assert(omg[n*nZ] < EPSILON);
		assert(omg[n*nZ+nZ-1] < EPSILON);
	}

	// Stuff for critical rayleigh check
	double tmpPrev[nN];
	double omgPrev[nN];
	double psiPrev[nN];
	for(int n=0; n<nN; ++n){
		tmpPrev[n] = tmp[32+n*nZ];
		psiPrev[n] = psi[32+n*nZ];
		omgPrev[n] = omg[32+n*nZ];
	}
	current = 0;
	int steps = 0;
	while (t<totalTime) {
		if(steps%500 == 0) {
			for(int n=1; n<nN; ++n){
				printf("%d: %e, %e, %e\n",n,
						std::log(std::abs(tmp[32+n*nZ])) - std::log(std::abs(tmpPrev[n])),
						std::log(std::abs(omg[32+n*nZ])) - std::log(std::abs(omgPrev[n])),
						std::log(std::abs(psi[32+n*nZ])) - std::log(std::abs(psiPrev[n])));
				tmpPrev[n] = tmp[32+n*nZ];
				psiPrev[n] = psi[32+n*nZ];
				omgPrev[n] = omg[32+n*nZ];
				/*
				printf("%d: %e, %e, %e\n", n, 
						tmp[nZ*n+32],
						omg[nZ*n+32],
						psi[nZ*n+32]);
						*/

			}
			/*
			printMaxOf(tmp, "tmp");
			printMaxOf(omg, "omg");
			printMaxOf(psi, "psi");
			//printMaxOf(dOmgdt+current*nN*nZ, "dOmgdt");
			std::printf(" \n");
			*/
		}
		steps++;
		computeLinearDerivatives(1);
		updateTmpAndOmg();
		solveForPsi();
		t+=dt;
		++current%=2;
	}	
}

int main() {
// Sim::Sim(int nZ, int nN, double dt, double Ra, double Pr, int a ,double timeBetweenSaves, bool modifydt, int current, double t, double totalTime
	Sim simulation = Sim(101, 51, 3.0e-6, 1e6, 0.5, 3, 1.5e-3, false, 0, 0, 3e-3);
	simulation.runNonLinear();
	//simulation.runLinear();

	// test_Sim_triDiagonalSolver();
	// test_Sim_dfdz2();
	// std::cout << test_Sim_dfdz() << std::endl;
}

