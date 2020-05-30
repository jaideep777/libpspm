//#include "solver.h"

#include <iostream>
#include <cmath>
#include <cassert>
#include <string>


using namespace std;


std::vector <double> seq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

// ~~~~~~~~~~~ SOLVER ~~~~~~~~~~~~~~~~~~~~~
template<class Model>
Solver<Model>::Solver(std::vector<double> xbreaks, PSPM_SolverType _method){
	xb = xbreaks[0];
	xm = xbreaks[xbreaks.size()-1];
	J  = xbreaks.size()-1;
	method = _method;

	x = xbreaks;

	X.resize(J);
	for (size_t i=0; i<J; ++i) X[i] = (xbreaks[i]+xbreaks[i+1])/2.0;
	
	h.resize(J);	// This will be used by FMU and updated by MMU, so pre-allocate
	for (size_t i=0; i<J; ++i) h[i] = xbreaks[i+1] - xbreaks[i];	
	
	switch (method){
		case SOLVER_FMU:	
			nx = J;						// X
			state_size = J;				// U
			state.resize(state_size,0);
			break;

		case SOLVER_MMU:
			nx = J;						// X
			state_size = J+1 + J;		// x, U
			state.resize(state_size,0);
			for (size_t i=0; i<J+1; ++i) state[i] = xbreaks[i];
			break;

		case SOLVER_CM:
			nx = J+1;					// x
			state_size = J+1 + J+1;		// x, u
			state.resize(state_size,0);
			for (size_t i=0; i<J+1; ++i) state[i] = xbreaks[i];
			break;

		case SOLVER_EBT:
			nx = J+1;					// [x0, xint]
			state_size = 1 + J + 1 + J;	// pi0, xint, N0, Nint
			state.resize(state_size,0);
			for (size_t i=0; i<J; ++i) state[1+i] = X[i];	// leave [0] for pi0 (= 0)
			X.insert(X.begin(), xb);	// x0 = xb + pi0/N0
			break;
	}


	rates.resize(state_size);
	
}

template<class Model>
void Solver<Model>::setModel(Model *M){
	mod = M;
}

	
template<class Model>
Solver<Model>::Solver(int _J, double _xb, double _xm, PSPM_SolverType _method) 
	: Solver(seq(_xb, _xm, _J+1), _method){
}


template<class Model>
const int Solver<Model>::size(){
	return state_size;
}

template<class Model>
const int Solver<Model>::xsize(){
	return nx;	
}

	//for (int i=0; i<J; ++i){
	//    x[i+1] = exp(log(0.01) + (i)*(log(xm)-log(0.01))/(J-1));
	//    h[i] = x[i+1]-x[i];
	//    X[i] = (x[i+1]+x[i])/2;
	//}



template<class Model>
const double * Solver<Model>::getX(){
	switch (method){
		case SOLVER_FMU:	return x.data(); 
		case SOLVER_MMU:	return x.data();
		case SOLVER_CM:		return x.data();
		case SOLVER_EBT:	return X.data();
	}
}


template<class Model>
void Solver<Model>::print(){
	string types[] = {"FMU", "MMU", "CM", "EBT"};
	std::cout << "Type: " << types[method] << std::endl;

	std::cout << "X (" << X.size() << "): ";
	for (int i=0; i<X.size(); ++i) cout << X[i] << " ";
	std::cout << std::endl;
	
	std::cout << "State (" << state.size() << "):\n";
	std::cout << "  X (" << state.size()-nx << "): ";
	for (int i=0; i<state.size()-nx; ++i) cout << state[i] << " ";
	std::cout << "\n  U (" << nx << "): ";
	for (int i=0; i<nx; ++i) cout << state[state.size()-nx + i] << " ";
	std::cout << std::endl;
}





template <class Model>
void Solver<Model>::initialize(){
	// state vector was initialized to 0 in Constrctor. Set non-zero elements here
	switch (method){
		case SOLVER_FMU:	
			for (size_t i=0; i<J; ++i)  state[i] = mod->calcIC(X[i]);
			break;

		case SOLVER_MMU:
			for (size_t i=0; i<J; ++i)  state[J+1 + i] = mod->calcIC(X[i]);
			//for (size_t i=0; i<J+1; ++i) uprev[i] = calcIC(x[i]);
			break;

		case SOLVER_CM:
			for (size_t i=0; i<J+1; ++i)  state[J+1 + i] = mod->calcIC(x[i]);
			break;

		case SOLVER_EBT:
			for (size_t i=0; i<J; ++i)  state[J+1 + 1+i] = mod->calcIC(X[1+i]);	// state[J+1+0]=0 (N0)
			break;
	}
}

template<class Model>
template<typename wFunc>
double Solver<Model>::integrate_x(wFunc w, int power){
	if (method == SOLVER_FMU || method == SOLVER_MMU){
		// integrate using midpoint quadrature rule
		double I=0;
		double * U = state.data();
		for (unsigned int i=0; i<X.size(); ++i){
			I += h[i]*w(X[i])*pow(U[i], power);  // TODO: Replace with std::transform after profiling
		}
		return I;
	}
	else{
		std::cout << "Only FMU and MMU are implemented\n";
		return 0;
	}
}


// ~~~~~~~~~~~ FMU Solver ~~~~~~~~~~~
double phi(double r){
	return max(max(0.0,min(2*r,1.0)),min(r,2.0));
}


template <class Model>
void Solver<Model>::calcRates_FMU(double t){

	double * U = state.data();

	vector <double> growthArray(J+1);
	for (int i=0; i<J+1; ++i) growthArray[i] = mod->growthRate(x[i], t, mod->evalEnv(x[i],t));

//	#define growth(i) growthRate(x[i], mod->evalEnv(x[i],t))
	#define growth(i) growthArray[i]

	// i=0
	double birthFlux = 0;
	for (int j=0; j<J; ++j) birthFlux += h[j]*mod->birthRate(X[j], t, mod->evalEnv(X[j],t))*U[j];
	
	vector <double> u(J+1);

	u[0] = birthFlux/(growth(0)+1e-12); // Q: is this correct? or g(X0,env)? 
	//cout << env.time << ": " << birthFlux/(growthRate(x[0],env)+1e-6) << endl;
	
	// i=1 (calc u1 assuming linear u(x) in first interval)
	u[1] = 2*U[0]-u[0];  // NOTE: for g(x) < 0 this can be calculated with upwind scheme 
	
	for (int i=2; i<J-1; ++i){ // dU[i] ~ u[i+1] <-- U[i],U[i-1], u[i] <-- U[i-1],U[i-2]
		if(growth(i) >=0){
			double rMinus = ((U[i]-U[i-1])/(x[i]-x[i-1]))/((U[i-1]-U[i-2]+1e-12)/(x[i-1]-x[i-2]));
			u[i] = U[i-1] + phi(rMinus)*(U[i-1]-U[i-2])*(x[i]-x[i-1])/(x[i+1]-x[i-1]); 
		}   
		else{
			double rPlus  = ((U[i]-U[i-1])/(x[i]-x[i-1]))/((U[i+1]-U[i]+1e-12)/(x[i+1]-x[i]));
			u[i] = U[i] - phi(rPlus)*(U[i+1]-U[i])*(x[i+1]-x[i])/(x[i+2]-x[i]); 
		}
	}
	
	u[J-1] = 2*U[J-2] - u[J-2];	// NOTE: for g(x) > 0 This can be calc with upwind scheme
	u[J] = 2*U[J-1] - u[J-1];

	for (int i=0; i<J; ++i){ // dU[i] ~ u[i+1] <-- U[i],U[i-1], u[i] <-- U[i-1],U[i-2]

		//if (i == 0) assert( fabs(growthRate(x[i],env)*u[i] - birthFlux) < 1e-11);
		rates[i] = -mod->mortalityRate(X[i], t, mod->evalEnv(X[i], t))*U[i] - (growth(i+1)*u[i+1] - growth(i)*u[i])/h[i];
		//f[i] = dU[i];
		
//		cout << growthRate(x[i],env)*u[i] << " " << growthRate(x[i-1],env)*u[i-1] << " = " << (growthRate(x[i],env)*u[i] - growthRate(x[i-1],env)*u[i-1]) << "\n--";
	}
//	cout << endl;	
//	cout << dU[0] << " " << dU[1] << endl;
	
}

	
//template<class Model> 
//void Solver<Model>::step_to(double tf){
	//double dt = 0.1;
	
	//for (int t=t_now; i<

//}
