#include <cassert>
#include <string>
#include "solver.h"
using namespace std;

void Solver::stepU_iCM(double t, vector<double> &S, vector<double> &dSdt, double dt){
		
	vector<double>::iterator its = S.begin() + n_statevars_system; // Skip system variables
	
	// 1. Take implicit step for U
	for (int s = 0; s<species_vec.size(); ++s){
		Species_Base* spp = species_vec[s];
		
		// [S S S x u x u x u x u x u a b c a b c a b c a b c a b c] <--- full SV for species is this
		//        ^ its, itr
		double *XU = &(*its); 
		int J = spp->J;

		std::vector<double> g_gx = spp->growthRateGradient(-1, spp->xb, t, env, control.cm_grad_dx);
		//std::cout << "g = " << g_gx[0] << ", gx = " << g_gx[1] << "\n";

		// Boundary u is not used in rate calcs per se, but needed in size integral. Hence update
		double gb = g_gx[0], growthGrad = g_gx[1];	
		double pe = spp->establishmentProbability(t, env);
		if (spp->birth_flux_in < 0){	
			double birthFlux = calcSpeciesBirthFlux(s,t) * pe;
			spp->set_ub(birthFlux/(gb+1e-20));
		}
		else{
			spp->calc_boundary_u(gb, pe); // this will set u in boundary cohort
		}

		// all cohorts
		for (int i=0; i<J; ++i){
			XU[2*i+0] += spp->growthRate(i, spp->getX(i), t, env)*dt;

			std::vector<double> g_gx = spp->growthRateGradient(i, spp->getX(i), t, env, control.cm_grad_dx);
			XU[2*i+1] /= 1 + g_gx[1]*dt + spp->mortalityRate(i, spp->getX(i), t, env)*dt;
		}
		
		its += J*(2+spp->n_extra_statevars);

	}
	
	assert(distance(S.begin(), its)==S.size());
}


