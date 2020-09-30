#include <algorithm>
#include <cassert>

#include "vector_insert.h"

template <class Model>
void Solver<Model>::calcRates_CM(double t, vector<double>&S, vector<double> &dSdt){
	auto is = createIterators_state(S);
	auto ir = createIterators_rates(dSdt);
	auto& itx = is.get("X");
	auto& itu = is.get("u");
	
	auto& itdx = ir.get("X");
	auto& itdu = ir.get("u");
	auto& itre = ir.get((varnames_extra.size()>0)? varnames_extra[0] : "X");	// dummy init
	
	for (is.begin(), ir.begin(); !is.end(); ++is, ++ir){
		double grad_dx = 0.001;
		
		double gxplus = mod->growthRate(*itx + grad_dx, t); 
		double gx     = mod->growthRate(*itx, t); 
		
		double growthGrad = (gxplus-gx)/grad_dx;

		*itdx =  gx;
		*itdu = -mod->mortalityRate(*itx, t)*(*itu) - growthGrad*(*itu);
		
		if (varnames_extra.size() > 0){
			auto it_returned = mod->calcRates_extra(t, *itx, itre);
			assert(distance(itre, it_returned) == varnames_extra.size());
		}
	}
	
	//double * x = &S[0];
	//double * u = &S[J+1];

	//double * dx = &dSdt[0];
	//double * du = &dSdt[J+1];

	//for (size_t i=0; i<J+1; ++i){
		
	//}

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
	cout << "xsize = " << xsize() << " " << J << endl;
	vector <int> at = {0, xsize()};
	vector <double> vals = {xb+3, -1};
	for (int i=0; i<varnames_extra.size(); ++i){
		at.push_back(2*xsize());
		vals.push_back(-2); // FIXME: get these from initStateExtra()
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
		state[J+1] = u0_in;
	}

}


template <class Model>
void Solver<Model>::removeCohort_CM(){
	// cohorts are x0, x1, x2, x3, ...xJ,  u0, u1, u2, u3, .... uJ 
	auto px = state.begin(); advance(px, 1); // point at x1
	auto pu = state.begin(); advance(pu, J+1+1); // point at u1
	auto last = state.begin(); advance(last, J-1+1); // point at xJ (1 past the last value to be considered)

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



