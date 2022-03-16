#include <cassert>
#include "solver.h"
using namespace std;

void Solver::calcRates_iFMU(double t, vector<double>::iterator S, vector<double>::iterator dSdt){
	vector<double>::iterator its = S    + n_statevars_system; // Skip system variables
	vector<double>::iterator itr = dSdt + n_statevars_system;
	
	for (int s = 0; s<species_vec.size(); ++s){
		auto spp = species_vec[s];
		
		its += spp->J;// skip u 
		for (int i=0; i<spp->J; ++i) *itr++ = 0; // set du/dt to 0 
	
		if (spp->n_extra_statevars > 0){
			auto itr_prev = itr;
			spp->getExtraRates(itr); // TODO/FIXME: Does calc of extra rates need t and env?
			assert(distance(itr_prev, itr) == spp->n_extra_statevars*spp->J);
			its += spp->n_extra_statevars*spp->J; 	
		}
	}

}

void Solver::stepU_iFMU(double t, vector<double> &S, vector<double> &dSdt, double dt){
		
	vector<double>::iterator its = S.begin() + n_statevars_system; // Skip system variables
	
	// 1. Take implicit step for U
	for (int s = 0; s<species_vec.size(); ++s){
		Species_Base* spp = species_vec[s];
		
		// [S S S u u u u u a b c a b c a b c a b c a b c] <--- full SV for species is this
		//        ^ its, itr
		double *U = &(*its); // Since FMU only has U in state, start of species is actually U
		int J = spp->J;			// xsize of species.
		vector<double> &h = spp->h;
		
		vector <double> growthArray(J);
		for (int i=0; i<J; ++i) growthArray[i] = spp->growthRate(i, spp->getX(i), t, env);
		
		double birthFlux;
		if (spp->birth_flux_in < 0){	
			birthFlux = calcSpeciesBirthFlux(s,t) * spp->establishmentProbability(t, env);
			//std::cout << "birthflux = " << calcSpeciesBirthFlux(s,t) << " * " << spp->establishmentProbability(t, env) << " = " << birthFlux << endl; 
		}
		else{
			birthFlux = spp->get_u0(t, env)*growthArray[0];
		}
		
		//cout << t << "\t" << birthFlux/growthArray[0] << " -> " << calcSpeciesBirthFlux(s,t)/growthArray[0] << "\n";
		double B0  = 1 + dt/h[0]*growthArray[0] + dt*spp->mortalityRate(0, spp->getX(0), t, env);
//		std::cout << "B0 = " << B0 << endl; 
		double C0 = spp->getU(0) + dt/h[0]*birthFlux;
		U[0] = C0/B0;
//		std::cout << "U0 = " << U[0] << endl; 

		for (int w = 1; w < J; ++w){
			double Aw = -growthArray[w-1]*dt/h[w];
			double Bw  = 1 + dt/h[w]*growthArray[w] + dt*spp->mortalityRate(w, spp->getX(w), t, env);
			double Cw = spp->getU(w);

			U[w] = (Cw - Aw*U[w-1])/Bw;
		}
		
		its += J*(1+spp->n_extra_statevars);

	}
	
}


