#include <cassert>

// ~~~~~~~~~~~ FMU Solver ~~~~~~~~~~~
double phi(double r){
	return max(max(0.0,min(2*r,1.0)),min(r,2.0));
}


template <class Model>
void Solver<Model>::calcRates_FMU(double t, vector<double> &U, vector<double> &dUdt){

	vector <double> growthArray(J+1);
	for (int i=0; i<J+1; ++i) growthArray[i] = mod->growthRate(x[i], t);

//	#define growth(i) growthRate(x[i], mod->evalEnv(x[i],t))
	#define growth(i) growthArray[i]

	vector <double> u(J+1);
	
	// i=0
	if (u0_in < 0){	
		double birthFlux = integrate_x([this](double z, double t){return mod->birthRate(z,t);}, t, U, 1);
		u[0] = birthFlux/(growth(0)+1e-12); // Q: is this correct? or g(X0,env)? 
	}
	else{
		u[0] = u0_in;
	}

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
		dUdt[i] = -mod->mortalityRate(X[i], t)*U[i] - (growth(i+1)*u[i+1] - growth(i)*u[i])/h[i];
	}
	
}

	

