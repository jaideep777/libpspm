#include <cassert>
#include "solver.h"
using namespace std;

// ~~~~~~~~~~~ FMU Solver ~~~~~~~~~~~
double phi(double r){
	return max(max(0.0,min(2*r,1.0)),min(r,2.0));
}


void Solver::calcRates_FMU(double t, vector<double>::iterator S, vector<double>::iterator dSdt){

	vector<double>::iterator its = S    + n_statevars_system; // Skip system variables
	vector<double>::iterator itr = dSdt + n_statevars_system;
	
	for (int s = 0; s<species_vec.size(); ++s){
		auto spp = species_vec[s];
		
		// [S S S u u u u u a b c a b c a b c a b c a b c] <--- full SV for species is this
		//        ^ its, itr
		double *U = &(*its); // Since FMU only has U in state, start of species is actually U
		double *dUdt = &(*itr); 
		int J = spp->J;			// xsize of species.
		vector<double> &x = spp->x;
		vector<double> &X = spp->X;
		vector<double> &h = spp->h;

		//spp->preComputeAllCohorts(t,env);
		
		vector <double> growthArray(J+1);
		growthArray[0] = spp->growthRate(-1, x[0], t, env);
		for (int i=1; i<J+1; ++i) growthArray[i] = spp->growthRateOffset(i-1, x[i], t, env);

	//	#define growth(i) growthRate(x[i], mod->evalEnv(x[i],t))
		//#define growth(i) growthArray[i]

		//for (int i=0; i<J+1; ++i){
			//std::cout << x[i] << " " << growthArray[i] << "\n";
		//}

		vector <double> u(J+1);
		
		// i=0
		if (spp->birth_flux_in < 0){	
			double birthFlux = calcSpeciesBirthFlux(s,t);
			//double birthFlux = integrate_x([this, s](int i, double t){
												//return species_vec[s]->birthRate(i,species_vec[s]->getX(i),t, env);
											//}, t, s);
											
			u[0] = birthFlux/(growthArray[0]+1e-12); // Q: is this correct? or g(X0,env)? - A: It is g(xb,env) - growthArray[] indexes x, so (0) is xb. Hence correct. 
		}
		else{
			u[0] = spp->get_u0(t, env);
			//double g = spp.mod->growthRate(spp.xb, current_time, env);
			//// --- debug ---
			//if (spp.bfin_is_u0in)
				//u[0] = spp.birth_flux_in;
			//else{
			//// -------------
				//double d = (g>0)? spp.birth_flux_in * spp.mod->establishmentProbability(current_time, env)/g  : 0;
				//u[0] = d; 
			//}
		}

		// i=1 (calc u1 assuming linear u(x) in first interval)
		u[1] = 2*U[0]-u[0];  // NOTE: for g(x) < 0 this can be calculated with upwind scheme 
		
		for (int i=2; i<J-1; ++i){ // dU[i] ~ u[i+1] <-- U[i],U[i-1], u[i] <-- U[i-1],U[i-2]
			if(growthArray[i] >=0){
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

		// [S S S u u u u u a b c a b c a b c a b c a b c] <--- full SV for species is this
		//        ^ itr, its. Therefore, advance 1 at a time while setting dUdt, and advance its by 3*5 after.
		for (int i=0; i<spp->J; ++i){ // dU[i] ~ u[i+1] <-- U[i],U[i-1], u[i] <-- U[i-1],U[i-2]
			dUdt[i] = -spp->mortalityRate(i,X[i], t, env)*U[i] - (growthArray[i+1]*u[i+1] - growthArray[i]*u[i])/h[i];
			++itr; ++its;
		}
		if (spp->n_extra_statevars > 0){
			auto itr_prev = itr;
			spp->getExtraRates(itr);
			assert(distance(itr_prev, itr) == spp->n_extra_statevars*spp->J);
			its += spp->n_extra_statevars*spp->J; 	
		}
			

	}
}

	

