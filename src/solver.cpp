#include "solver.h"

#include <iostream>
#include <cmath>
#include <cassert>
#include <string>


using namespace std;

//double Environment::eval(double x){
//    return value;
//}


//double growthRate(double x, Environment &env){
//    double E = env.eval(x);
//    double t = env.time;
//    double a = 0.16+0.22*exp(-0.225*t*t);
//    return 0.225*(1-x*x)*(E/(1+E*E))*t*(1+a*a)/a;
//}

//double mortalityRate(double x, Environment &env){
//    double E = env.eval(x);
//    double t = env.time;
//    double a = 0.16+0.22*exp(-0.225*t*t);
//    return 1.35*t*E/a;
//}

//double birthRate(double x, Environment &env){
//    double E = env.eval(x);
//    double t = env.time;
//    double oneplusa = 1.16+0.22*exp(-0.225*t*t);
//    double a = 0.16+0.22*exp(-0.225*t*t);
//    double n1 = 0.225*t*x*x*(1-x)*(1-x)*E/(1+E)/(1+E)*oneplusa*oneplusa/a;
//    double n2 = (1+exp(-0.225*t*t))/(61-88*log(2)+(38*log(2)-79.0/3)*exp(-0.225*t*t));
//    return n1*n2;
//}


//double computeEnvironment(vector <double>& x, vector <double>& u, vector<double>& h){
//    vector <double> w(x.size());
//    for (unsigned int i=0; i<x.size(); ++i){
//        if (x[i] <= 1.0/3){
//            w[i] = 1;		
//        }
//        else if (x[i] > 1.0/3 && x[i] <= 2.0/3){
//            w[i] = pow(2-3*x[i], 3)*(54*x[i]*x[i]-27*x[i]+4);
//        }
//        else {
//            w[i]= 0;
//        }
//    } 
	
//    // integrate using midpoint quadrature rule
//    double I=0;
//    for (unsigned int i=0; i<x.size(); ++i) I += h[i]*w[i]*u[i];
	
//    return I;
	
//}
	
//vector<double> initialize(vector <double>&x){
//    vector<double> u(x.size());
//    for (unsigned int i=0; i<x.size(); ++i){
//        u[i] = pow(1-x[i],2)/pow(1+x[i],4) + (1-x[i])/pow(1+x[i],3);
//    }
//    return u;
//}


std::vector <double> seq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

// ~~~~~~~~~~~ SOLVER ~~~~~~~~~~~~~~~~~~~~~
Solver::Solver(std::vector<double> xbreaks, PSPM_SolverType _method){
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


Solver::Solver(int _J, double _xb, double _xm, PSPM_SolverType _method) 
	: Solver(seq(_xb, _xm, _J+1), _method){
}


const int Solver::size(){
	return state_size;
}

const int Solver::xsize(){
	return nx;	
}

	//for (int i=0; i<J; ++i){
	//    x[i+1] = exp(log(0.01) + (i)*(log(xm)-log(0.01))/(J-1));
	//    h[i] = x[i+1]-x[i];
	//    X[i] = (x[i+1]+x[i])/2;
	//}



const double * Solver::getX(){
	switch (method){
		case SOLVER_FMU:	return x.data();
		case SOLVER_MMU:	return x.data();
		case SOLVER_CM:		return x.data();
		case SOLVER_EBT:	return X.data();
	}
}

//void Solver::initialize(std::vector <double>& u0){
	//// state vector was initialized to 0 in Constrctor. Set non-zero elements here
	//switch (method){
		//case SOLVER_FMU:	
			//for (size_t i=0; i<J; ++i)  state[i] = u0[i];
			//break;

		//case SOLVER_MMU:
			//for (size_t i=0; i<J; ++i)  state[J+1 + i] = u0[i];
			//break;

		//case SOLVER_CM:
			//for (size_t i=0; i<J; ++i)  state[J+1 + i] = u0[i];
			//break;

		//case SOLVER_EBT:
			//for (size_t i=0; i<J; ++i)  state[J+1 + 1+i] = u0[i];	// leave [J+1+0] for N0
			//break;
	//}
//}




void Solver::print(){
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



//double phi(double r){
//    return max(max(0.0,min(2*r,1.0)),min(r,2.0));
//}

//void Solver::setState(const double * y){
//    for (int i=0; i<J; ++i) U[i] = y[i];
//}

//void Solver::getRates(Environment &env, double * f){
	
//    vector <double> growthArray(J+1);
//    for (int i=0; i<J+1; ++i) growthArray[i] = growthRate(x[i],env);

////	#define growth(i) growthRate(x[i],env)
//    #define growth(i) growthArray[i]

//    // i=0
//    double birthFlux = 0;
//    for (int j=0; j<J; ++j) birthFlux += h[j]*birthRate(X[j],env)*U[j];
//    u[0] = birthFlux/(growth(0)+1e-12); // Q: is this correct? or g(X0,env)? 
//    //cout << env.time << ": " << birthFlux/(growthRate(x[0],env)+1e-6) << endl;
	
//    // i=1 (calc u1 assuming linear u(x) in first interval)
//    u[1] = 2*U[0]-u[0]; 
	
//    for (int i=1; i<J; ++i){ // dU[i] ~ u[i+1] <-- U[i],U[i-1], u[i] <-- U[i-1],U[i-2]
//        if (i >= 2){
//            if(growth(i) >=0){
//                double rMinus = (U[i]-U[i-1])*(x[i]-x[i-2])/((U[i-1]-U[i-2])*(x[i-1]-x[i-2]));
//                u[i] = U[i-1] + phi(rMinus)*(U[i-1]-U[i-2])*(x[i]-x[i-1])/(x[i+1]-x[i-1]); 
//            }   
//            else{
//                double rPlus  = (U[i-1]-U[i])*(x[i]-x[i+2])/((U[i]-U[i+1])*(x[i-1]-x[i]));
//                u[i] = U[i] - phi(rPlus)*(U[i+1]-U[i])*(x[i+1]-x[i])/(x[i+2]-x[i]); 
//            }
//        }
//    }
	
//    u[J] = 2*U[J-1] - u[J-1];

//    for (int i=0; i<J; ++i){ // dU[i] ~ u[i+1] <-- U[i],U[i-1], u[i] <-- U[i-1],U[i-2]

//        if (i == 0) assert( fabs(growthRate(x[i],env)*u[i] - birthFlux) < 1e-11);
//        dU[i] = -mortalityRate(X[i], env)*U[i] - (growth(i+1)*u[i+1] - growth(i)*u[i])/h[i];
//        f[i] = dU[i];
		
////		cout << growthRate(x[i],env)*u[i] << " " << growthRate(x[i-1],env)*u[i-1] << " = " << (growthRate(x[i],env)*u[i] - growthRate(x[i-1],env)*u[i-1]) << "\n--";
//    }
////	cout << endl;	
////	cout << dU[0] << " " << dU[1] << endl;
	
//}




