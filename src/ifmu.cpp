#include <cassert>
#include "solver.h"

void Solver::step_iFMU(double t, vector<double> &S, double dt){
	
	vector<double>::iterator its = S.begin()    + n_statevars_system; // Skip system variables
	
	for (int s = 0; s<species_vec.size(); ++s){
		auto spp = species_vec[s];
		
		// [S S S u u u u u a b c a b c a b c a b c a b c] <--- full SV for species is this
		//        ^ its, itr
		double *U = &(*its); // Since FMU only has U in state, start of species is actually U
		int J = spp->J;			// xsize of species.
		vector<double> &h = spp->h;
		
		vector <double> growthArray(J);
		for (int i=0; i<J; ++i) growthArray[i] = spp->growthRate(i, spp->getX(i), t, env);
		
		double birthFlux;
		if (spp->birth_flux_in < 0){	
			birthFlux = calcSpeciesBirthFlux(s,t);
		}
		else{
			birthFlux = spp->get_u0(t, env)*growthArray[0];
		}
		
		cout << t << "\t" << birthFlux/growthArray[0] << " -> " << calcSpeciesBirthFlux(s,t)/growthArray[0] << "\n";
		double B0  = 1 + dt/h[0]*growthArray[0] + dt*spp->mortalityRate(0, spp->getX(0), t, env);
		double C0 = spp->getU(0) + dt/h[0]*birthFlux;
		U[0] = C0/B0;

		for (int w = 1; w < J; ++w){
			double Aw = -growthArray[w-1]*dt/h[w];
			double Bw  = 1 + dt/h[w]*growthArray[w] + dt*spp->mortalityRate(w, spp->getX(w), t, env);
			double Cw = spp->getU(w);

			U[w] = (Cw - Aw*U[w-1])/Bw;
		}
		
		for (int i=0; i<J; ++i) spp->setU(i, U[i]); // copy new U to cohorts
	}
	
}



