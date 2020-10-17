#include <algorithm>
#include <cassert>

#include "vector_insert.h"

template<class Model, class Environment>
void Solver<Model,Environment>::calcRates_CM(double t, vector<double>&S, vector<double> &dSdt){

	for (auto& spp : species_vec){	
		auto is = spp.get_iterators(S);
		auto ir = spp.get_iterators(dSdt);
		auto& itx = is.get("X");
		auto& itu = is.get("u");
		auto& itse = is.get((spp.varnames_extra.size()>0)? spp.varnames_extra[0] : "X");	// dummy init //FIXME: commented line below does not work. Explore why
		//auto& itse = (varnames_extra.size()>0)? is.get(varnames_extra[0]) : S.begin(); 
		
		auto& itdx = ir.get("X");
		auto& itdu = ir.get("u");
		auto& itre = ir.get((spp.varnames_extra.size()>0)? spp.varnames_extra[0] : "X");	// dummy init
		//auto& itre = (varnames_extra.size()>0)? ir.get(varnames_extra[0]) : dSdt.begin(); 
		
		for (is.begin(), ir.begin(); !is.end(); ++is, ++ir){
			double grad_dx = 1e-6;
			
			double gxplus = spp.mod->growthRate(*itx + grad_dx, t, env); 
			double gx     = spp.mod->growthRate(*itx, t, env); 
			
			double growthGrad = (gxplus-gx)/grad_dx;

			*itdx =  gx;
			*itdu = -(spp.mod->mortalityRate(*itx, t, env) + growthGrad); //*(*itu);
			
			if (spp.varnames_extra.size() > 0){
				auto it_returned = spp.mod->calcRates_extra(t, *itx, itse, itre);
				assert(distance(itre, it_returned) == spp.varnames_extra.size());
			}
		}
	}
}


template<class Model, class Environment>
double Solver<Model,Environment>::calc_u0_CM(){
	//// function to iterate
	//auto f = [this](double utry){
	//    // set u0 to given (trial) value
	//    state[xsize()] = utry;
	//    // recompute environment based on new state
	//    mod->computeEnv(current_time, state, this);
	//    // calculate birthflux by trapezoidal integration (under new environment)
	//    double birthFlux = integrate_x([this](double z, double t){return mod->birthRate(z,t);}, current_time, state, 1);
		
	//    double unext = birthFlux/mod->growthRate(xb, current_time);
	//    return unext;
	//};

	//double u0 = state[xsize()+1]; // initialize with u0 = u1
	//// iterate
	//double err = 100;
	//while(err > 1e-6){
	//    double u1 = f(u0);
	//    err = abs(u1 - u0);
	//    u0 = u1;
	//}
	//state[xsize()] = u0;
	////cout << "u0 = " << u0 << endl;	
	//return state[xsize()];
	return 0;
}


template<class Model, class Environment>
void Solver<Model,Environment>::addCohort_CM(){
	// Create index and value vectors to insert into state
	vector <int> at;
	vector <double> vals;
	int n_insertions = 0;	// keep track of cumulative number of insertions so that species start_index can be updated

	for (auto& spp : species_vec){
		int spp_start_index_new = spp.start_index + n_insertions; // calculate new start_index for species and store it

		at.push_back(spp.start_index + 0);		// TODO: should these be generalized to use iteratorset?
		at.push_back(spp.start_index + spp.J);
		vals.push_back(spp.xb);
		vals.push_back(-100);
		n_insertions += 2;
		vector <double> ex = spp.mod->initStateExtra(spp.xb, current_time);
		for (int i=0; i<spp.varnames_extra.size(); ++i){
			at.push_back(spp.start_index + 2*spp.J);
			vals.push_back(ex[i]);
			++n_insertions;	
		}
		
		spp.start_index = spp_start_index_new;	// update start_index of species
		++spp.J;								// increment J to reflect new system size
		setupLayout(spp);						// rebuild parameters for IteratorSet
	}

	vector_insert(state, at, vals);		// insert new cohorts into state
	rates.resize(state.size(), -999);	// resize rates to correct size

	// initialize the newly inserted cohorts
	for (auto&spp : species_vec){
		// set u0 to correct value
		if (spp.birth_flux_in < 0){
			calc_u0_CM(); // this internally sets u0 in state
		}
		else {
			double g = spp.mod->growthRate(spp.xb, current_time, env);
			// --- debug ---
			if (spp.bfin_is_u0in)
				state[spp.start_index + spp.J] = log(spp.birth_flux_in);
			else
			// -------------
				state[spp.start_index + spp.J] = (g>0)? log(spp.birth_flux_in * spp.mod->establishmentProbability(current_time)/g)  :  log(0); 
		}
	}
}


template<class Model, class Environment>
void Solver<Model,Environment>::removeCohort_CM(){
	// TODO: this should be made truely generic - maybe after making of indexset
	//// cohorts are x0, x1, ...xJ-1,  u0, u1, .... uJ-1, a0,b0,c0,..., a1,b1,c1,...,aJ-1,bJ-1,cJ-1
	//auto px = state.begin(); advance(px, 1); // point at x1
	//auto pu = state.begin(); advance(pu, xsize()+1); // point at u1
	//auto last = state.begin(); advance(last, xsize()-1); // point at xJ (1 past the last value to be considered)

	////cout << *px << " " << *last << " " << *pu << endl;
	vector<double> flags(state.size(), 0.0);

	int n_removals = 0;
	for (auto& spp : species_vec){
		int spp_start_index_new = spp.start_index - n_removals;
		auto iset = spp.get_iterators(state);
		auto iset_f = spp.get_iterators(flags);
		
		auto& itx = iset.get("X");
		++iset; ++iset_f; // point to second element
		auto iset_remove = iset_f;
		double dx_min = *std::next(itx) - *std::prev(itx);
		for (; iset.dist < iset.size-1; ++iset, ++iset_f){	// go till second-last element
			double dx = *std::next(itx) - *std::prev(itx);
			if (dx < dx_min){
				dx_min = dx;
				iset_remove = iset_f;
			}
			//cout << dx << " " << dx_min << " " << *remove_x << " " << *remove_u << endl;
		}
	
		// mark flags to remove
		auto iset_vec_f = iset_remove.get();
		for (int i=0; i<iset_vec_f.size(); ++i){
			*iset_vec_f[i] = 1.0;
			++n_removals;
		}	
	
		spp.start_index = spp_start_index_new;
		--spp.J; // reduce species size
		setupLayout(spp);	
	}
	
	// remove all flagged elements from state
	auto pred = [this, &flags](const double& s) -> bool {
		return flags[&s - (const double*)&this->state[0]] == 1.0; 
	};

	auto pend = std::remove_if(state.begin(), state.end(), pred);
	state.erase(pend, state.end());
	////for (auto z : state) cout << z << " "; cout << endl;

}



