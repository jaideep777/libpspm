#include <cassert>
#include <string>
#include "solver.h"
using namespace std;

// void Solver::calcRates_iEBT(double t, vector<double>::iterator S, vector<double>::iterator dSdt){
// 	vector<double>::iterator its = S    + n_statevars_system; // Skip system variables
// 	vector<double>::iterator itr = dSdt + n_statevars_system;
	
// 	for (int s = 0; s<species_vec.size(); ++s){
// 		auto spp = species_vec[s];
		
// 		its += (n_statevars_internal)*spp->J; // skip x and u 
// 		for (int i=0; i<n_statevars_internal*spp->J; ++i) *itr++ = 0; // set dx/dt and du/dt to 0 
	
// 		if (spp->n_extra_statevars > 0){
// 			auto itr_prev = itr;
// 			spp->getExtraRates(itr); // TODO/FIXME: Does calc of extra rates need t and env?
// 			assert(distance(itr_prev, itr) == spp->n_extra_statevars*spp->J);
// 			its += spp->n_extra_statevars*spp->J; 	
// 		}
// 	}

// }


void Solver::stepU_iEBT(double t, vector<double> &S, vector<double> &dSdt, double dt){
		
	vector<double>::iterator its = S.begin() + n_statevars_system; // Skip system variables
	
	// 1. Take implicit step for U
	for (int s = 0; s<species_vec.size(); ++s){
		Species_Base* spp = species_vec[s];
		
		// [S S S x u x u x u x u x u a b c a b c a b c a b c a b c] <--- full SV for species is this
		//        ^ its, itr
		double *XU = &(*its); 
		int J = spp->J;

		double   pi0  =  spp->getX(spp->J-1);	 // last cohort is pi0, N0
		double   N0   =  spp->getU(spp->J-1);
		//std::cout << "pi = " << pi0 << ", N0 = " << N0 << "\n";

		std::vector<double> g_gx = spp->growthRateGradient(-1, spp->xb, t, env, control.ebt_grad_dx);
		std::vector<double> m_mx = spp->mortalityRateGradient(-1, spp->xb, t, env, control.ebt_grad_dx);
		//std::cout << "g = " << g_gx[0] << ", gx = " << g_gx[1] << "\n";
		//std::cout << "m = " << m_mx[0] << ", mx = " << m_mx[1] << "\n";

		double mb = m_mx[0], mortGrad = m_mx[1], gb = g_gx[0], growthGrad = g_gx[1];	
		
		double birthFlux;
		double pe = spp->establishmentProbability(t, env);
		if (spp->birth_flux_in < 0){	
			birthFlux = calcSpeciesBirthFlux(s,t) * pe;
		}
		else{
			double u0 = spp->calc_boundary_u(gb, pe);
			birthFlux = u0*gb;
		}

		// internal cohorts
		for (int i=0; i<J-1; ++i){
			XU[2*i+0] += spp->growthRate(i, spp->getX(i), t, env)*dt;
			XU[2*i+1] /= 1+spp->mortalityRate(i, spp->getX(i), t, env)*dt;
		}

		// pi0 cohort
		double a1 = 1 + mb*dt;
		double b1 = mortGrad*dt;
		double c1 = birthFlux*dt + N0;
		double a2 = -gb*dt;
		double b2 = 1 - growthGrad*dt + mb*dt;
		double c2 = pi0;

		double pinew = (a2*c1-a1*c2)/(a2*b1-a1*b2);
		double unew = (b2*c1-b1*c2)/(b2*a1-b1*a2);

		if (pinew < 0) throw std::runtime_error("pi0 < 0: "+std::to_string(pinew));
		if (unew < 0) throw std::runtime_error("u0 < 0: "+std::to_string(unew));
		
		XU[2*(J-1)+0] = pinew;
		XU[2*(J-1)+1] = unew;

		its += J*(2+spp->n_extra_statevars);

	}
	
	assert(distance(S.begin(), its)==S.size());
}


