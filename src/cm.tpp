#include <algorithm>
#include <cassert>

template <class Model>
void Solver<Model>::calcRates_CM(double t, vector<double>&S, vector<double> &dSdt){
	double * x = &S[0];
	double * u = &S[J+1];

	double * dx = &dSdt[0];
	double * du = &dSdt[J+1];

	for (size_t i=0; i<J+1; ++i){
		dx[i] =  mod->growthRate(x[i], t);
		
		double grad_dx = 0.001;
		double growthGrad = (mod->growthRate(x[i]+0.001, t) - mod->growthRate(x[i], t))/0.001;
		
		du[i] = -mod->mortalityRate(x[i], t)*u[i] - growthGrad*u[i];
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
	auto p_x  = state.begin();
	auto p_u  = state.begin(); advance(p_u, J+1); 

	state.insert(p_u, -1); // insert new u0 (dummy value) BEFORE p_u. (Now both iterators are invalid)
	state.insert(state.begin(), xb); // this inserts xb at the 1st position	
	++J;	// increment J to reflect new system size

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



