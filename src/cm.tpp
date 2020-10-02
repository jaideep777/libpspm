#include <algorithm>
#include <cassert>

#include "vector_insert.h"

template <class Model>
void Solver<Model>::calcRates_CM(double t, vector<double>&S, vector<double> &dSdt){
	auto is = createIterators_state(S);
	auto ir = createIterators_rates(dSdt);
	auto& itx = is.get("X");
	auto& itu = is.get("u");
	auto& itse = is.get((varnames_extra.size()>0)? varnames_extra[0] : "X");	// dummy init
	
	auto& itdx = ir.get("X");
	auto& itdu = ir.get("u");
	auto& itre = ir.get((varnames_extra.size()>0)? varnames_extra[0] : "X");	// dummy init
	
	for (is.begin(), ir.begin(); !is.end(); ++is, ++ir){
		double grad_dx = 1e-6;
		
		double gxplus = mod->growthRate(*itx + grad_dx, t); 
		double gx     = mod->growthRate(*itx, t); 
		
		double growthGrad = (gxplus-gx)/grad_dx;

		*itdx =  gx;
		*itdu = -(mod->mortalityRate(*itx, t) + growthGrad); //*(*itu);
		
		if (varnames_extra.size() > 0){
			auto it_returned = mod->calcRates_extra(t, *itx, itse, itre);
			assert(distance(itre, it_returned) == varnames_extra.size());
		}
	}
	
}


template <class Model>
double Solver<Model>::calc_u0_CM(){
	// function to iterate
	auto f = [this](double utry){
		// set u0 to given (trial) value
		state[J+1] = utry;
		// recompute environment based on new state
		mod->computeEnv(current_time, state, this);
		// calculate birthflux by trapezoidal integration (under new environment)
		double birthFlux = integrate_x([this](double z, double t){return mod->birthRate(z,t);}, current_time, state, 1);
		
		double unext = birthFlux/mod->growthRate(xb, current_time);
		return unext;
	};

	double u0 = state[J+2]; // initialize with u0 = u1
	// iterate
	double err = 100;
	while(err > 1e-6){
		double u1 = f(u0);
		err = abs(u1 - u0);
		u0 = u1;
	}
	state[J+1] = u0;
	//cout << "u0 = " << u0 << endl;	
	return state[J+1];
}


template <class Model>
void Solver<Model>::addCohort_CM(){
	//auto p_x  = state.begin();
	//auto p_u  = state.begin(); advance(p_u, J+1); 
	//cout << "xsize = " << xsize() << " " << J << endl;
	vector <int> at = {0, xsize()};
	vector <double> vals = {xb, -1};
	vector <double> ex = mod->initStateExtra(xb, current_time);
	for (int i=0; i<varnames_extra.size(); ++i){
		at.push_back(2*xsize());
		vals.push_back(ex[i]); 
	}
	
	vector_insert(state, at, vals);
	rates.resize(state.size(), -999);
	//state.insert(p_u, -1); // insert new u0 (dummy value) BEFORE p_u. (Now both iterators are invalid)
	//state.insert(state.begin(), xb); // this inserts xb at the 1st position	
	++J;	// increment J to reflect new system size
	setupLayout();  // rebuild parameters for IteratorSet

	// set u0 to correct value
	if (u0_in < 0){
		calc_u0_CM(); // this internally sets u0 in state
	}
	else {
		double g = mod->growthRate(xb, current_time);
		//cout << "g = " << g << "\n";
		state[J+1] = (g>0)? log(u0_in*mod->establishmentProbability(current_time)/g)  :  log(0); //FIXME: set to 0 if g()<0
	}

}


template <class Model>
void Solver<Model>::removeCohort_CM(){
	// cohorts are x0, x1, x2, x3, ...xJ,  u0, u1, u2, u3, .... uJ 
	auto px = state.begin(); advance(px, 1); // point at x1
	auto pu = state.begin(); advance(pu, xsize()+1); // point at u1
	auto last = state.begin(); advance(last, xsize()-1); // point at xJ (1 past the last value to be considered)

	//cout << *px << " " << *last << " " << *pu << endl;

	double dx_min = *std::next(px) - *std::prev(px);
	auto remove_x = px;
	auto remove_u = pu;
	while(px != last){
		double dx = *std::next(px) - *std::prev(px);
		if (dx < dx_min){
			dx_min = dx;
			remove_x = px;
			remove_u = pu;
		}
		++px; ++pu;
		//cout << dx << " " << dx_min << " " << *remove_x << " " << *remove_u << endl;
	}
	state.erase(remove_u);    // remove farther element first so remove_x remains valid
	state.erase(remove_x);
	--J;

	//for (auto z : state) cout << z << " "; cout << endl;

}



