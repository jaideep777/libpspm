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
Solver<Model>::Solver(std::vector<double> xbreaks, PSPM_SolverType _method) : odeStepper(0){
	xb = xbreaks[0];
	xm = xbreaks[xbreaks.size()-1];
	J  = xbreaks.size()-1;
	method = _method;

	x = xbreaks;

	X.resize(J);
	for (size_t i=0; i<J; ++i) X[i] = (xbreaks[i]+xbreaks[i+1])/2.0;
	
	h.resize(J);	// This will be used by FMU and updated by MMU, EBT, so pre-allocate
	for (size_t i=0; i<J; ++i) h[i] = xbreaks[i+1] - xbreaks[i];	
	
	if (method == SOLVER_FMU){	
		//nx = J;						// X
		state_size = J;				// U
		state.resize(state_size,0);
	}

	if (method == SOLVER_MMU){
		//nx = J;						// X
		state_size = J+1 + J;		// x, U
		state.resize(state_size,0);
		for (size_t i=0; i<J+1; ++i) state[i] = xbreaks[i];
	}

	if (method == SOLVER_CM){
		//nx = J+1;					// x
		state_size = J+1 + J+1;		// x, u
		state.resize(state_size,0);
		for (size_t i=0; i<J+1; ++i) state[i] = xbreaks[i];
	}

	if (method == SOLVER_EBT){
		//nx = J+1;					// [x0, xint]
		state_size = 1 + J + 1 + J;	// pi0, xint, N0, Nint
		state.resize(state_size,0);
		for (size_t i=0; i<J; ++i) state[1+i] = X[i];	// leave [0] for pi0 (= 0)
		X.insert(X.begin(), xb);	// x0 = xb + pi0/N0
	}


	rates.resize(state_size);

}

template<class Model>
void Solver<Model>::setModel(Model *M){
	mod = M;
}


template<class Model>
void Solver<Model>::setInputNewbornDensity(double input_u0){
	u0_in = input_u0;
}


template<class Model>
Solver<Model>::Solver(int _J, double _xb, double _xm, PSPM_SolverType _method) 
	: Solver(seq(_xb, _xm, _J+1), _method){
}


template<class Model>
const int Solver<Model>::size(){
	return state.size();
}

template<class Model>
const int Solver<Model>::xsize(){
	if (method == SOLVER_FMU) return J;	
	if (method == SOLVER_MMU) return J;  
	if (method == SOLVER_CM ) return J+1;
	if (method == SOLVER_EBT) return J+1;
}

	//for (int i=0; i<J; ++i){
	//    x[i+1] = exp(log(0.01) + (i)*(log(xm)-log(0.01))/(J-1));
	//    h[i] = x[i+1]-x[i];
	//    X[i] = (x[i+1]+x[i])/2;
	//}



template<class Model>
const double * Solver<Model>::getX(){
	if (method == SOLVER_FMU)	return x.data();  
	if (method == SOLVER_MMU)	return x.data();
	if (method == SOLVER_CM )	return x.data();
	if (method == SOLVER_EBT)	return X.data();
}

template<class Model>
vector<double> Solver<Model>::getx(){
	return x;
}


template<class Model>
void Solver<Model>::print(){
	string types[] = {"FMU", "MMU", "CM", "EBT"};
	std::cout << "Type: " << types[method] << std::endl;

	std::cout << "X (" << X.size() << "): ";
	for (int i=0; i<X.size(); ++i) cout << X[i] << " ";
	std::cout << std::endl;
	
	std::cout << "State (" << state.size() << "):\n";
	std::cout << "  X (" << state.size()-xsize() << "): ";
	for (int i=0; i<state.size()-xsize(); ++i) cout << state[i] << " ";
	std::cout << "\n  U (" << xsize() << "): ";
	for (int i=0; i<xsize(); ++i) cout << state[state.size()-xsize() + i] << " ";
	std::cout << std::endl;
}





template <class Model>
void Solver<Model>::initialize(){
	// state vector was initialized to 0 in Constrctor. Set non-zero elements here
	if (method == SOLVER_FMU){
		for (size_t i=0; i<J; ++i)  state[i] = mod->calcIC(X[i]);
	}
	if (method == SOLVER_MMU){
		for (size_t i=0; i<J; ++i)  state[J+1 + i] = mod->calcIC(X[i]);
		//for (size_t i=0; i<J+1; ++i) uprev[i] = calcIC(x[i]);
	}
	if (method == SOLVER_CM){
		for (size_t i=0; i<J+1; ++i)  state[J+1 + i] = mod->calcIC(x[i]);
	}
	if (method == SOLVER_EBT){
		for (size_t i=0; i<J; ++i)  state[J+1 + 1+i] = mod->calcIC(X[1+i])*h[i];	// state[J+1+0]=0 (N0)
	}
}

template<class Model>
template<typename wFunc>
double Solver<Model>::integrate_x(wFunc w, double t, vector<double>&S, int power){
	//cout << " | " <<  t << " " << mod->evalEnv(0,t) << " ";
	if (method == SOLVER_FMU || method == SOLVER_MMU){
		// integrate using midpoint quadrature rule
		double I=0;
		double * U = S.data();
		for (unsigned int i=0; i<X.size(); ++i){
			I += h[i]*w(X[i], t)*pow(U[i], power);  // TODO: Replace with std::transform after profiling
		}
		return I;
	}
	else if (method == SOLVER_EBT){
		// integrate using EBT rule (sum over cohorts)
		double   pi0  =  S[0];
		double * xint = &S[1];
		double   N0   =  S[J+1];
		double * Nint = &S[J+2];
		
		double x0 = xb + pi0/(N0+1e-12); 
		
		double I = w(x0, t)*N0;
		for (int i=0; i<J; ++i) I += w(xint[i], t)*Nint[i];
		
		return I;
	}
	else if (method == SOLVER_CM){
		// integrate using trapezoidal rule TODO: Modify to avoid double computation of w(x)
		double * px = &S[0];
		double * pu = &S[J+1];
		double I = 0;
		for (int i=0; i<J; ++i){
			I += (px[i+1]-px[i])*(w(px[i+1], t)*pu[i+1]+w(px[i], t)*pu[i]);
		}
		return I*0.5;
	}
	else{
		std::cout << "Only FMU and MMU are implemented\n";
		return 0;
	}
}


// current_time is updated by the ODE solver at every (internal) step
template<class Model>
void Solver<Model>::step_to(double tstop){
	// do nothing if tstop is <= current_time
	if (tstop <= current_time) return;
	
	if (method == SOLVER_FMU){	
		auto derivs = [this](double t, vector<double> &S, vector<double> &dSdt){
			mod->computeEnv(t, S, this);
			this->calcRates_FMU(t, S, dSdt);
		};
		
		odeStepper.Step_to(tstop, current_time, state, derivs); // state = [U]
	}
	if (method == SOLVER_MMU){
	}
	if (method == SOLVER_EBT){
		auto derivs = [this](double t, vector<double> &S, vector<double> &dSdt){
			mod->computeEnv(t, S, this);
			this->calcRates_EBT(t, S, dSdt);
		};
		
		// integrate 
		odeStepper.Step_to(tstop, current_time, state, derivs); // state = [pi0, Xint, N0, Nint]
		
		// update cohorts
		removeDeadCohorts_EBT();
		if (state[J+1] > 0) addCohort_EBT();  // Add new cohort if N0 > 0. Add after removing dead ones otherwise this will also be removed. 
		
		// update variables based on new state
		// X, x, h, etc  ==> Maybe not necessary

	}
	if (method == SOLVER_CM){
		auto derivs = [this](double t, vector<double> &S, vector<double> &dSdt){
			mod->computeEnv(t, S, this);
			this->calcRates_CM(t, S, dSdt);
		};
		
		// integrate 
		odeStepper.Step_to(tstop, current_time, state, derivs); // state = [pi0, Xint, N0, Nint]

		// update cohorts
		addCohort_CM();		// add before so that it becomes boundary cohort and first internal cohort can be (potentially) removed
		removeCohort_CM();

	}
}


template<class Model>
double Solver<Model>::newborns_out(){
	// update Environment from latest state
	mod->computeEnv(current_time, state, this);
	// calculate birthflux 
	double birthFlux = integrate_x([this](double z, double t){return mod->birthRate(z,t);}, current_time, state, 1);
	return birthFlux;
}


template<class Model>
double Solver<Model>::u0_out(){
	return newborns_out()/mod->growthRate(xb, current_time);
}


#include "mu.tpp"
#include "ebt.tpp"
#include "cm.tpp"


