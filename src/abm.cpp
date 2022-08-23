#include <cassert>
#include "solver.h"
using namespace std;

void Solver::stepABM(double t, double dt){
		
	// 1. Take implicit step for U
	for (int s = 0; s<species_vec.size(); ++s){
		Species_Base* spp = species_vec[s];
		
		vector<double> growthArray(spp->J);
		for (int i=0; i<spp->J; ++i){
			// Calc growth rates
			growthArray[i] = spp->growthRate(i, spp->getX(i), t, env);
		}

		// implement fecundity
		double pe = spp->establishmentProbability(t, env);
		double gb = spp->growthRate(-1, spp->xb, t, env);
		double birthFlux;
		if (spp->birth_flux_in < 0){	
			birthFlux = calcSpeciesBirthFlux(s,t) * pe;
		}
		else{
			double ub = spp->get_boundary_u(); // backup boundary cohort's density because it will be overwritten by calc
			double u0 = spp->calc_boundary_u(gb, pe);
			birthFlux = u0*gb;
			spp->set_ub(ub);
		}

		double noff = birthFlux*dt/spp->get_boundary_u();  // number of offspring = birthflux / density of each superindividual
		//cout << "noff: " << spp->noff_abm << " | " << noff << " " << birthFlux << " " << dt << " " << spp->get_boundary_u() << endl;
		spp->noff_abm += noff;
		
		// implement mortality
		for (int i=0; i<spp->J; ++i){
			double mortRate = spp->mortalityRate(i, spp->getX(i), t, env);
			double survival_prob = exp(-mortRate*dt);
			//cout << "survival prob = " << survival_prob << "\n";
			if (double(rand())/RAND_MAX >= survival_prob) spp->markCohortForRemoval(i);
		}
		spp->removeMarkedCohorts();

		// implement growth
		for (int i=0; i<spp->J; ++i){
			spp->setX(i, spp->getX(i) + growthArray[i]*dt);
		}


		// step extra variables
		// these will be in order abc abc abc... 
		// extra rates
		vector<double>::iterator its = state.begin() + n_statevars_system; // Skip system variables
		vector<double>::iterator itr = rates.begin() + n_statevars_system;
		if (spp->n_extra_statevars > 0){
			auto itr_prev = itr;
			auto its_prev = its;

			spp->copyCohortsExtraToState(its); // this will increment its
			spp->getExtraRates(itr);      // this will increment itr
			assert(distance(itr_prev, itr) == spp->n_extra_statevars*spp->J);
			assert(distance(its_prev, its) == spp->n_extra_statevars*spp->J);

			itr = itr_prev; its = its_prev; // so bring them back
			
			for (int i=0; i< spp->n_extra_statevars*spp->J; ++i){
				*its += (*itr)*dt;
				++its; ++itr;
			}
			itr = itr_prev; its = its_prev; // so bring them back

			spp->copyExtraStateToCohorts(its);
		}
		
		// add recruits once they have accumulated to > 1
		if (spp->noff_abm > 1){
			int nadd = int(spp->noff_abm);	
			spp->addCohort(nadd);
			spp->noff_abm -= nadd;
		}
		// also add a recruit in anticipation if the whole populaiton is dead
		if (spp->J == 0){
			spp->addCohort(1);
			//--spp->noff_abm;
		}
		
		// resize state - this matters only if extra istate is specified
		resizeStateFromSpecies();
	}
	
	//print();
	
}


