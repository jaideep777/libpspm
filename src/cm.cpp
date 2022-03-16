#include <algorithm>
#include <cassert>

#include "solver.h"
//#include "vector_insert.h"
using namespace std;

// state must be copied to cohorts before calling this function
void Solver::calcRates_CM(double t, vector<double>::iterator S, vector<double>::iterator dSdt){

	//auto ss = species_vec[0];	
	//cout << "svec: " << t << " " << S[0] << " " << S[1] << " " << S[2] << " " << S[3] << "\n";
	//cout << "coho: " << t << " " << ss->getX(0) << " " << ss->getU(0) << " " << ss->getX(1) << " " << ss->getU(1) << "\n";
	
	vector<double>::iterator its = S    + n_statevars_system; // Skip system variables
	vector<double>::iterator itr = dSdt + n_statevars_system;

	for (int s = 0; s<species_vec.size(); ++s){
		Species_Base* spp = species_vec[s];
		//cout << "calcRates(species = " << s << ")\n";
		for (int i=0; i<spp->J; ++i){
			double x = spp->getX(i);
			vector<double> g_gx = spp->growthRateGradient(i, x, t, env, control.cm_grad_dx);  // FIXME: x can go
			double dx =  g_gx[0];
			double du = -(spp->mortalityRate(i, x, t, env) + g_gx[1]); 
			if (!use_log_densities) du *= spp->getU(i);
		
			*itr++ = dx;
			*itr++ = du;	
			
			its+=2;
		}
		if (spp->n_extra_statevars > 0){
			auto itr_prev = itr;
			spp->getExtraRates(itr);
			assert(distance(itr_prev, itr) == spp->n_extra_statevars*spp->J); 
			its += spp->n_extra_statevars*spp->J; 	
		}
		//cout << "---\n";
	}
}


void Solver::addCohort_CM(){

	for (auto& spp : species_vec){
		spp->get_u0(current_time, env);  // init density of boundary cohort
		spp->initBoundaryCohort(current_time, env); // init extra state variables and birth time of the boundary cohort
		spp->addCohort();	// introduce copy of boundary cohort in system
	}
	
	resizeStateFromSpecies();
	copyCohortsToState();
}


void Solver::removeCohort_CM(){

	for (auto& spp : species_vec){
		if (spp->J > control.max_cohorts) // remove a cohort if number of cohorts in the species exceeds a threshold
			spp->removeDensestCohort();
	}
	
	resizeStateFromSpecies();
	copyCohortsToState();
	
}



//template<class Model, class Environment>
//double Solver<Model,Environment>::calc_u0_CM(){
	////// function to iterate
	////auto f = [this](double utry){
	////    // set u0 to given (trial) value
	////    state[xsize()] = utry;
	////    // recompute environment based on new state
	////    mod->computeEnv(current_time, state, this);
	////    // calculate birthflux by trapezoidal integration (under new environment)
	////    double birthFlux = integrate_x([this](double z, double t){return mod->birthRate(z,t);}, current_time, state, 1);
		
	////    double unext = birthFlux/mod->growthRate(xb, current_time);
	////    return unext;
	////};

	////double u0 = state[xsize()+1]; // initialize with u0 = u1
	////// iterate
	////double err = 100;
	////while(err > 1e-6){
	////    double u1 = f(u0);
	////    err = abs(u1 - u0);
	////    u0 = u1;
	////}
	////state[xsize()] = u0;
	//////cout << "u0 = " << u0 << endl;	
	////return state[xsize()];
	//return 0;
//}

