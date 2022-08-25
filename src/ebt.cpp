#include <algorithm>
#include <cassert>

#include <solver.h>
using namespace std;


void Solver::calcRates_EBT(double t, vector<double>::iterator S, vector<double>::iterator dSdt){

	vector<double>::iterator its = S    + n_statevars_system; // Skip system variables
	vector<double>::iterator itr = dSdt + n_statevars_system;

	for (int s = 0; s < species_vec.size(); ++s){
		Species_Base * spp = species_vec[s];

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

		for (int i=0; i<spp->J-1; ++i){	// go down to the second last cohort (exclude boundary cohort)
			double dx = spp->growthRate(i, spp->getX(i), t, env);							// dx/dt
			double du = -spp->mortalityRate(i, spp->getX(i), t, env) * spp->getU(i);		// du/dt
			//std::cout << "S/C = " << s << "/" << i << " " << spp->getX(i) << " " << dx << " " << du << "\n";
			*itr++ = dx;
			*itr++ = du;
			its += 2;
		}

		// dpi0/dt and dN0/dt
		//std::cout << "S/C = " << s << "/" << "b" << " | pi0/N0 = " << pi0 << " " << N0;
		//if (pi0 <= 0) pi0 = 1e-40;
		//if (N0  <= 0) N0  = 1e-40;
		double dN0  = -mb*N0 - mortGrad*pi0 + birthFlux;
		double dpi0 = gb*N0 + growthGrad*pi0 - mb*pi0;
		//std::cout << " | dpi0/dN0 = " << dpi0 << " " << dN0 << " | mx/gx/mb/gb = " << m_mx[1] << " " << g_gx[1] << " " << m_mx[0] << " " << g_gx[0] << "\n";
		*itr++ = dpi0;
		*itr++ = dN0;
		its += 2;

		if (spp->n_extra_statevars > 0){
			auto itr_prev = itr;
			spp->getExtraRates(itr);
			assert(distance(itr_prev, itr) == spp->n_extra_statevars*spp->J); 
			its += spp->n_extra_statevars*spp->J; 	
		}
		
	}

}



void Solver::addCohort_EBT(){
	// Q: is updateEnv needed here (and in cm.cpp) to init cumulative vars?
	// Maybe not, see logic in abm.cpp
	for (auto spp : species_vec){
		// 1. internalize the pi0-cohort (this cohort's birth time would have been already set when it was inserted)
		// - get pi0, N0 from last cohort
		double   pi0  =  spp->getX(spp->J-1);
		double   N0   =  spp->getU(spp->J-1);

		// - update the recently internalized pi0-cohort with actual x0 value
		double x0 = spp->xb + pi0/(N0+1e-12);
		spp->setX(spp->J-1, x0);
		spp->setU(spp->J-1, N0);

		// 2. insert a new cohort (copy of boundary cohort, to be the new pi0-cohort)
		spp->initBoundaryCohort(current_time, env);	// update initial extra state and birth-time of boundary cohort
		spp->addCohort(); // introduce copy of boundary cohort into species

		// 3. set x,u of the new cohort to 0,0, thus marking it as the new pi0-cohort
		spp->setX(spp->J-1, 0); 
		spp->setU(spp->J-1, 0);
	}

	// 4. reset state from cohorts
	resizeStateFromSpecies(); 
	copyCohortsToState();
	
}


void Solver::removeDeadCohorts_EBT(){

	for (auto spp : species_vec){
		spp->removeDeadCohorts(control.ebt_ucut);	
	}

	resizeStateFromSpecies();
	copyCohortsToState();

}


