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

	int s = 0;
	for (auto spp : species_vec){	
		cout << "calcRates(species = " << s << ")\n";
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
		++s;
		cout << "---\n";
	}
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


void Solver::addCohort_CM(){
	
	//copyStateToCohorts(state.begin());
	//int state_size_new = n_statevars_system;
	for (auto& spp : species_vec){
		spp->get_u0(current_time, env);  // init density of boundary cohort
		spp->initBoundaryCohort(current_time, env); // init extra state variables and birth time of the boundary cohort
		spp->addCohort();	// introduce copy of boundary cohort in system
		//state_size_new += spp->J*(n_statevars_internal + spp->n_extra_statevars);
	}
	
	//state.resize(state_size_new, -999);
	//rates.resize(state_size_new, -999);
	resizeStateFromSpecies();
	copyCohortsToState();

	// if (birth_flux_in < 0) calc_u0_cm...

	// old bad code...
	//int spp_start_index_new = spp.start_index + n_insertions; // calculate new start_index for species and store it

		//at.push_back(spp.start_index + 0);		// TODO: should these be generalized to use iteratorset?
		//at.push_back(spp.start_index + spp.J);
		//vals.push_back(spp.xb);
		//vals.push_back(-100);
		//n_insertions += 2;
		//vector <double> ex = spp.mod->initStateExtra(spp.xb, current_time, env);
		//for (int i=0; i<spp.varnames_extra.size(); ++i){
			//at.push_back(spp.start_index + 2*spp.J);
			//vals.push_back(ex[i]);
			//++n_insertions;	
		//}
		
		//spp.start_index = spp_start_index_new;	// update start_index of species
		//++spp.J;								// increment J to reflect new system size
		//setupLayout(spp);						// rebuild parameters for IteratorSet
	//}

	//vector_insert(state, at, vals);		// insert new cohorts into state
	//rates.resize(state.size(), -999);	// resize rates to correct size

	//// initialize the newly inserted cohorts
	//for (auto&spp : species_vec){
	//    // set u0 to correct value
	//    if (spp.birth_flux_in < 0){
	//        calc_u0_CM(); // this internally sets u0 in state
	//    }
	//    else {
	//        double g = spp.mod->growthRate(spp.xb, current_time, env);
	//        // --- debug ---
	//        if (spp.bfin_is_u0in)
	//            state[spp.start_index + spp.J] = (use_log_densities)? log(spp.birth_flux_in) : spp.birth_flux_in;
	//        else{
	//        // -------------
	//            double d = (g>0)? spp.birth_flux_in * spp.mod->establishmentProbability(current_time, env)/g  : 0;
	//            state[spp.start_index + spp.J] = (use_log_densities)? log(d) : d; 
	//        }
	//    }
	//}
}


//template<class Model, class Environment>
void Solver::removeCohort_CM(){
	// TODO: this should be made truely generic - maybe after making of indexset
	//// cohorts are x0, x1, ...xJ-1,  u0, u1, .... uJ-1, a0,b0,c0,..., a1,b1,c1,...,aJ-1,bJ-1,cJ-1
	//auto px = state.begin(); advance(px, 1); // point at x1
	//auto pu = state.begin(); advance(pu, xsize()+1); // point at u1
	//auto last = state.begin(); advance(last, xsize()-1); // point at xJ (1 past the last value to be considered)

	////cout << *px << " " << *last << " " << *pu << endl;
	//vector<double> flags(state.size(), 0.0);

	//copyStateToCohorts(state.begin());
	//int n_removals = 0;
	for (auto& spp : species_vec){
		if (spp->J > control.max_cohorts) // remove a cohort if number of cohorts in the species exceeds a threshold
			spp->removeDensestCohort();
	}
	
	resizeStateFromSpecies();
	copyCohortsToState();
	
	//int spp_start_index_new = spp.start_index - n_removals;
		//auto iset = spp.get_iterators(state);
		//auto iset_f = spp.get_iterators(flags);
		
		//auto& itx = iset.get("X");
		//++iset; ++iset_f; // point to second element
		//auto iset_remove = iset_f;
		//double dx_min = *std::next(itx) - *std::prev(itx);
		//for (; iset.dist < iset.size-1; ++iset, ++iset_f){	// go till second-last element
			//double dx = *std::next(itx) - *std::prev(itx);
			//if (dx < dx_min){
				//dx_min = dx;
				//iset_remove = iset_f;
			//}
			////cout << dx << " " << dx_min << " " << *remove_x << " " << *remove_u << endl;
		//}
	
		//// mark flags to remove
		//auto iset_vec_f = iset_remove.get();
		//for (int i=0; i<iset_vec_f.size(); ++i){
			//*iset_vec_f[i] = 1.0;
			//++n_removals;
		//}	
	
		//spp.start_index = spp_start_index_new;
		//--spp.J; // reduce species size
		//setupLayout(spp);	
	//}
	
	//// remove all flagged elements from state
	//auto pred = [this, &flags](const double& s) -> bool {
		//return flags[&s - (const double*)&this->state[0]] == 1.0; 
	//};

	//auto pend = std::remove_if(state.begin(), state.end(), pred);
	//state.erase(pend, state.end());
	//////for (auto z : state) cout << z << " "; cout << endl;

}


//template<class Model, class Environment>
//void Solver<Model,Environment>::removeDenseCohorts_CM(){
	//// TODO: this should be made truely generic - maybe after making of indexset
	////// cohorts are x0, x1, ...xJ-1,  u0, u1, .... uJ-1, a0,b0,c0,..., a1,b1,c1,...,aJ-1,bJ-1,cJ-1
	////auto px = state.begin(); advance(px, 1); // point at x1
	////auto pu = state.begin(); advance(pu, xsize()+1); // point at u1
	////auto last = state.begin(); advance(last, xsize()-1); // point at xJ (1 past the last value to be considered)

	//////cout << *px << " " << *last << " " << *pu << endl;
	//vector<double> flags(state.size(), 0.0);

	//int n_removals = 0;
	//for (auto& spp : species_vec){
		////if (spp.J < 500) continue;   // remove a cohort if number of cohorts in the species exceeds a threshold

		//int spp_start_index_new = spp.start_index - n_removals;  // start_index must be shifted left by number of elements removed so far
		//auto iset = spp.get_iterators(state);
		//auto iset_f = spp.get_iterators(flags);
		
		//auto& itx = iset.get("X");
		
		//iset.begin(); iset_f.begin();
		//++iset; ++iset_f; // point to second element (skip boundary cohort)

		////int n_cohorts_removed = 0;
		//for (; iset.dist < iset.size-1; ++++iset, ++++iset_f){	// go till second-last element
			//double dist = std::min(abs(*std::next(itx) - *itx), abs(*itx - *std::prev(itx)));
			//if (dist < 1e-4){
				//--spp.J; //++n_cohorts_removed;
				//// mark flags to remove
				//auto iset_vec_f = iset_f.get();  // vector of iterators to remove. FIXME: How to ensure that extra variables dont come up in this vector which are not part of state?
				//for (int i=0; i<iset_vec_f.size(); ++i){
					//*iset_vec_f[i] = 1.0;
					//++n_removals;
				//}
			//}
		//}
	
		//spp.start_index = spp_start_index_new;
		////spp.J -= n_cohorts_removed; // reduce species size
		//setupLayout(spp);	
	//}
	
	//// remove all flagged elements from state
	//auto pred = [this, &flags](const double& s) -> bool {
		//return flags[&s - (const double*)&this->state[0]] == 1.0; 
	//};

	//auto pend = std::remove_if(state.begin(), state.end(), pred);
	//state.erase(pend, state.end());
	//////for (auto z : state) cout << z << " "; cout << endl;

//}


