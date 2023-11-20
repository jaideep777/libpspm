#include <cassert>
#include "solver.h"
using namespace std;

// ~~~~~~~~~~~ FMU Solver ~~~~~~~~~~~
double phi(double r){
	return max(max(0.0,min(2*r,1.0)),min(r,2.0));
}


/// This solver uses the upwind scheme to discretize the McKendrick von-Foerster PDE.
/// \f[
/// \frac{\partial u}{\partial t} = - \frac{\partial}{\partial x}(gu) - \mu u
/// \f]
/// This yields:
/// The grid is set up from the vector of edges \f$x\f$, as specified during initialization.
/// Each point in the grid corresponds to a cell edge 
/// \f[
/// 
/// \f]
void Solver::calcRates_FMU(double t, vector<double>::iterator S, vector<double>::iterator dSdt){

	vector<double>::iterator its = S    + n_statevars_system; // Skip system variables
	vector<double>::iterator itr = dSdt + n_statevars_system;
	
	for (int s = 0; s<species_vec.size(); ++s){
		Species_Base* spp = species_vec[s];
		
		// [S S S u u u u u a b c a b c a b c a b c a b c] <--- full SV for species is this
		//        ^ its, itr
		double *U = &(*its); // Since FMU only has U in state, start of species is actually U
		double *dUdt = &(*itr); 
		int J = spp->J;			// xsize of species.
		vector<vector<double>> &x = spp->x;
		vector<vector<double>> &X = spp->X;
		vector<vector<double>> &h = spp->h;

		vector <double> growthArray(J+1);
		growthArray[0] = spp->growthRate(-1, x[0], t, env); // growth rate of boundary cohort
		for (int i=1; i<J+1; ++i) growthArray[i] = spp->growthRateOffset(i-1, x[i], t, env);

		vector <double> u(J+1);
		
		// i=0 -- set u0 from boundary condition
		double pe = spp->establishmentProbability(t, env);
		if (spp->birth_flux_in < 0){	
			double birthFlux = calcSpeciesBirthFlux(s,t) * pe;
			u[0] = birthFlux/(growthArray[0]+1e-12); // Q: is this correct? or g(X0,env)? - A: It is g(xb,env) - growthArray[] indexes x, so (0) is xb. Hence correct. 
		}
		else{
			u[0] = spp->calc_boundary_u(growthArray[0], pe);
		}

		// i=1 (calc u1 assuming linear u(x) in first interval) // NO: Horrible idea - leads to negative densities when u0 is very high
		// calc u[1] using first order method: u[1]=U[0] OR 2nd order method: u[1] = 2*U[0]-u[0]
		double u1_plus = fmax(2*U[0]-u[0], 0); // U[0] for constant interpolation, 2*U[0]-u[0] for linear
		u[1] = (growthArray[0] >= 0)? u1_plus : U[1]; 
		
		// i = 2 -- J-2
		for (int i=2; i<J-1; ++i){ // dU[i] ~ u[i+1] <-- U[i],U[i-1], u[i] <-- U[i-1],U[i-2]
			if(growthArray[i] >=0){
				//double rMinus = ((U[i]-U[i-1])/(x[i]-x[i-1]))/((U[i-1]-U[i-2]+1e-12)/(x[i-1]-x[i-2]));
				//u[i] = U[i-1] + phi(rMinus)*(U[i-1]-U[i-2])*(x[i]-x[i-1])/(x[i+1]-x[i-1]); 
				double r_down = (U[i]-U[i-1]) / (x[i+1]-x[i-1]); // factor of 2 not included because it will be cancelled in the ratio
				double r_up   = (U[i-1]-U[i-2]) / (x[i]-x[i-2]);
				double r = r_down/r_up;
				u[i] = U[i-1] + phi(r)*(U[i-1]-U[i-2])*(x[i]-x[i-1])/(x[i]-x[i-2]); 
			}   
			else{
				//throw std::runtime_error("");
				//double rPlus  = ((U[i]-U[i-1])/(x[i]-x[i-1]))/((U[i+1]-U[i]+1e-12)/(x[i+1]-x[i]));
				//u[i] = U[i] - phi(rPlus)*(U[i+1]-U[i])*(x[i+1]-x[i])/(x[i+2]-x[i]); 
				double r_down = (U[i]-U[i-1]) / (x[i+1]-x[i-1]); // factor of 2 not included because it will be cancelled in the ratio 
				double r_up   = (U[i+1]-U[i]) / (x[i+2]-x[i]);
				double r = r_down/r_up;
				u[i] = U[i] - phi(r)*(U[i+1]-U[i])*(x[i+1]-x[i])/(x[i+2]-x[i]); 
			}
		}
		
		// i = J-1 and i = J
		u[J-1] = (growthArray[J-1]>=0)? U[J-2] : U[J-1]; //2*U[J-2] - u[J-2];	// NOTE: for g(x) > 0 This can be calc with upwind scheme
		u[J] = U[J-1]; //2*U[J-1] - u[J-1];

		// core rates
		// [S S S u u u u u a b c a b c a b c a b c a b c] <--- full SV for species is this
		//        ^ itr, its. Therefore, advance 1 at a time while setting dUdt, and advance its by 3*5 after.
		for (int i=0; i<spp->J; ++i){ // dU[i] ~ u[i+1] <-- U[i],U[i-1], u[i] <-- U[i-1],U[i-2]
			dUdt[i] = -spp->mortalityRate(i,X[i], t, env)*U[i] - (growthArray[i+1]*u[i+1] - growthArray[i]*u[i])/h[i];
			++itr; ++its;
		}
		
		// extra rates
		if (spp->n_extra_statevars > 0){
			auto itr_prev = itr;
			spp->getExtraRates(itr);
			assert(distance(itr_prev, itr) == spp->n_extra_statevars*spp->J);
			its += spp->n_extra_statevars*spp->J; 	
		}
			

	}
}

	

