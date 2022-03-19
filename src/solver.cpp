#include "solver.h"
#include "cubic_spline.h"

#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include <algorithm>
#include <functional>


using namespace std;


inline std::vector <double> seq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

inline std::vector <double> logseq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = exp(log(from) + i*(log(to)-log(from))/(len-1));
	return x;
}

inline std::vector <double> diff(vector <double> breaks){
	std::vector<double> mids(breaks.size()-1);
	for (size_t i=0; i<mids.size(); ++i) mids[i] = (breaks[i]+breaks[i+1])/2;
	return mids;
}


// ~~~~~~~~~~~ SOLVER ~~~~~~~~~~~~~~~~~~~~~


Solver::Solver(PSPM_SolverType _method, string ode_method) : odeStepper(ode_method, 0, 1e-6, 1e-6) {
	method = _method;
}


void Solver::addSpecies(std::vector<double> xbreaks, Species_Base* s, int n_extra_vars, double input_birth_flux){
	s->set_inputBirthFlux(input_birth_flux);
	s->n_extra_statevars = n_extra_vars;

	if (method == SOLVER_FMU || method == SOLVER_IFMU){
		std::sort(xbreaks.begin(), xbreaks.end(), std::less<double>());  // sort cohorts ascending for FMU
		s->xb = *xbreaks.begin();
	}
	else {
		std::sort(xbreaks.begin(), xbreaks.end(), std::greater<double>());  // sort cohorts descending for CM, EBT
		s->xb = *xbreaks.rbegin();
	}

	int J;
	if      (method == SOLVER_FMU)  J = xbreaks.size()-1;	
	else if (method == SOLVER_IFMU) J = xbreaks.size()-1;	
	else if (method == SOLVER_MMU)  J = xbreaks.size()-1;  
	else if (method == SOLVER_CM )  J = xbreaks.size();
	else if (method == SOLVER_EBT)  J = xbreaks.size();
	else    throw std::runtime_error("Unsupported method");

	s->x = xbreaks;
	s->resize(J);

	// FMU has only 1 internal state variable (x), reset have 2 (x,u)
	bool cond = (method == SOLVER_FMU || method == SOLVER_IFMU);
	n_statevars_internal = (cond)? 1:2;
	
	int state_size = s->J*(n_statevars_internal + n_extra_vars);    // x/u and extra vars

	//s->start_index = state.size();	// New species will be appended to the end of state vector
	state.resize(state.size()+state_size);  // This will resize state for all species additions, but this in only initialization so its okay.
	rates.resize(rates.size()+state_size);

	species_vec.push_back(s);

	if (method == SOLVER_FMU){	
		s->X.resize(s->J);
		for (size_t i=0; i<s->J; ++i) s->X[i] = (s->x[i]+s->x[i+1])/2.0;
		
		s->h.resize(s->J);	// This will be used only by FMU
		for (size_t i=0; i<s->J; ++i) s->h[i] = s->x[i+1] - s->x[i];	
	}

	if (method == SOLVER_IFMU){
		s->X.resize(s->J);
		s->h.resize(s->J);	// This will be used only by FMU
		
		if (control.ifmu_centered_grids){ // grids are labelled by the grid center	
			for (size_t i=0; i<s->J; ++i) s->X[i] = (s->x[i]+s->x[i+1])/2.0;
			for (size_t i=0; i<s->J; ++i) s->h[i] = s->x[i+1] - s->x[i];	
			// for centered grids, also set xb to X[0] rather than x[0]. 
			s->xb = s->X[0];
		}
		else{  // grids are labelled by the lower edge (this ensures that recruitment happens exactly at xb)
			for (size_t i=0; i<s->J; ++i) s->X[i] = s->x[i];
			for (size_t i=0; i<s->J; ++i) s->h[i] = s->x[i+1] - s->x[i];	
			s->xb = s->X[0];
		}
	}

}


void Solver::addSpecies(int _J, double _xb, double _xm, bool log_breaks, Species_Base* s,
							   int n_extra_vars, double input_birth_flux){
	vector<double> xbreaks;
	if (log_breaks) xbreaks = logseq(_xb, _xm, _J+1);
	else            xbreaks =    seq(_xb, _xm, _J+1);

	addSpecies(xbreaks, s, n_extra_vars, input_birth_flux);
}


void Solver::addSystemVariables(int _s){
	n_statevars_system = _s;
	state.resize(state.size() + _s);
	rates.resize(state.size() + _s);
}


void Solver::resetState(double t0){  // FIXME: This is currently redundant, and needs to be improved with reset of both state and cohorts for a true reset of state
	current_time = t0;
	odeStepper.reset(t0, control.ode_eps, control.ode_eps); // = RKCK45<vector<double>> (0, control.ode_eps, control.ode_initial_step_size);  // this is a cheap operation, but this will empty the internal containers, which will then be (automatically) resized at next 1st ODE step. Maybe add a reset function to the ODE stepper? 

	std::fill(state.begin(), state.end(), 0); 
	std::fill(rates.begin(), rates.end(), -999); // DEBUG
}


void Solver::resizeStateFromSpecies(){
	int state_size_new = n_statevars_system;
	for (auto& spp : species_vec){
		state_size_new += spp->J*(n_statevars_internal + spp->n_extra_statevars);
	}
	
	state.resize(state_size_new, -999);
	rates.resize(state_size_new, -999);
}


void Solver::setEnvironment(EnvironmentBase * _env){
	env = _env;
}


////Species<Model>* Solver<Model,Environment>::get_species(int id){
//    return &species_vec[id];
//}


int Solver::n_species(){
	return species_vec.size();
}


double Solver::maxSize(){
	double maxsize = 0;
	for (auto spp : species_vec) maxsize = std::max(maxsize, spp->get_maxSize());
	return maxsize;
}




void Solver::print(){
	std::cout << ">> SOLVER \n";
	string types[] = {"FMU", "MMU", "CM", "EBT", "Implicit FMU"};
	std::cout << "+ Type: " << types[method] << std::endl;

	std::cout << "+ State size = " << state.size() << "\n";
	std::cout << "+ Rates size = " << rates.size() << "\n";
	std::cout << "+ Species:\n";
	std::cout.flush();
	for (int i=0; i<species_vec.size(); ++i) {
		std::cout << "Sp (" << i << "):\n";
		species_vec[i]->print();
	}
	
}


// Layout of the state vector is as follows:
//  ------------------------------------------------------------
// | x : u | x : u | x : u | a : b : c | a : b : c | a : b : c |
//  ------------------------------------------------------------
//  In above layout, the internal variables x and u are tightly
//  packed first, followed by user variables a, b, c. 
//  This arrangement is cache friendly because:
//  When setting rates, typically a,b,c are calculated one 
//  after the other for each x.
// TODO: should this take t0 as an argument, instead of setting to 0? 
void Solver::initialize(){
	
	vector<double>::iterator it = state.begin() + n_statevars_system; // TODO: replace with init_sState() 
	
	for (int k=0; k<species_vec.size(); ++k){
		Species_Base* s = species_vec[k];
	
		// Boundary cohort is not in state, but used as a reference.	
		s->set_xb(s->xb); // set x of boundary cohort 
		s->set_ub(0);     // set initial density of boundary cohort to 0.
		
		// set birth time for each cohort to current_time
		for (int i=0; i<s->J; ++i) s->set_birthTime(i, current_time);

		// set x, u for all cohorts
		if (method == SOLVER_FMU || method == SOLVER_IFMU){
			for (size_t i=0; i<s->J; ++i){
				double X = s->X[i];			// for FMU, X is the midpoint of the cell edges
				double U = s->init_density(i, X, env); 
				s->setX(i,X); 
				s->setU(i,U);
				*it++ = U;		// u in state (only)
			}
		}
		if (method == SOLVER_CM){
			for (size_t i=0; i<s->J; ++i){
				double X = s->x[i];
				double U = s->init_density(i, X, env); 
				s->setX(i,X); 
				s->setU(i,U);
				*it++ = X;									// x in state
				*it++ = (use_log_densities)? log(U) : U;	// u in state 
			}
		}
		if (method == SOLVER_EBT){
			// x, u for internal cohorts in state and it cohorts
			for (size_t i=0; i<s->J-1; ++i){
				double X = (s->x[i]+s->x[i+1])/2.0;			
				double U = s->init_density(i, X, env)*(s->x[i]-s->x[i+1]); 
				s->setX(i,X); 
				s->setU(i,U);
				*it++ = X;	// x in state
				*it++ = U;	// u in state 
			}
			// set pi0, N0 as x, u for the last cohort. This scheme allows using this last cohort with xb+pi0 in integrals etc 
			*it++ = 0; *it++ = 0;
			s->setX(s->J-1,0); // TODO: should this be set to xb for init_state and set to 0 again later? maybe not, as init_state is not expected to be x dependent
			s->setU(s->J-1,0); 		
		}

		// initialize extra state for each cohort and copy it to state
		s->initAndCopyExtraState(current_time, env, it);
		//if (s->n_extra_statevars > 0){  // FIXME: maybe redundant
			//auto it_prev = it;
			//s->init_ExtraState(it);  // this also inits the extra state of boundary cohort, but without advancing the iterator
			//assert(distance(it_prev, it) == s->n_extra_statevars*s->J); 
		//}

	}
}


void Solver::copyStateToCohorts(std::vector<double>::iterator state_begin){
	if (debug) std::cout << "state ---> cohorts\n";
	std::vector<double>::iterator it = state_begin + n_statevars_system; // no need to copy system state
	
	for (int k=0; k<species_vec.size(); ++k){
		Species_Base* s = species_vec[k];
		
		s->set_xb(s->xb); // Important: "touch" boundary cohort to trigger a precompute. Note that setSize() triggers it.
		if (method == SOLVER_FMU || method == SOLVER_IFMU){
			for (size_t i=0; i<s->J; ++i){
				s->setU(i,*it++);
			}
		}
		if (method == SOLVER_CM){
			for (size_t i=0; i<s->J; ++i){
				double X = *it++;	// get x from state
				double U = *it++;	// get u from state
			    U = (use_log_densities)? exp(U) : U;	// u in cohorts 
				s->setX(i,X); 
				s->setU(i,U);
			}
		}
		if (method == SOLVER_EBT){
			// x, u for boundary and internal cohorts
			for (size_t i=0; i<s->J; ++i){
				double X = *it++; 
				double U = *it++;
				s->setX(i,X); 
				s->setU(i,U);
			}
		}

		if (s->n_extra_statevars > 0){  // If extra state variables have been requested, initialize them
			auto it_prev = it;
			s->copyExtraStateToCohorts(it);
			assert(distance(it_prev, it) == s->n_extra_statevars*s->J); 
		}
	}
}


void Solver::copyCohortsToState(){
	if (debug) std::cout << "state <--- cohorts\n";
	vector<double>::iterator it = state.begin() + n_statevars_system; // no need to copy system state
	
	for (int k=0; k<species_vec.size(); ++k){
		Species_Base* s = species_vec[k];
		
		if (method == SOLVER_FMU || method == SOLVER_IFMU){
			for (size_t i=0; i<s->J; ++i){
				*it++ = s->getU(i);
			}
		}
		if (method == SOLVER_CM){
			for (size_t i=0; i<s->J; ++i){
				double X = s->getX(i); 
				double U = s->getU(i);
			    U = (use_log_densities)? log(U) : U;	// log(u) in state 
				*it++ = X;	// set x to state
				*it++ = U;	// set u to state
			}
		}
		if (method == SOLVER_EBT){
			// x, u for boundary and internal cohorts
			for (size_t i=0; i<s->J; ++i){
				double X = s->getX(i); 
				double U = s->getU(i);
				*it++ = X; 
				*it++ = U;
			}
		}

		if (s->n_extra_statevars > 0){  // If extra state variables have been requested, initialize them
			auto it_prev = it;
			s->copyCohortsExtraToState(it);
			assert(distance(it_prev, it) == s->n_extra_statevars*s->J); 
		}
	}
}



void Solver::step_to(double tstop){
	auto func = [](double t){};
	step_to(tstop, func);
}


void Solver::updateEnv(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt){
	if (debug) std::cout << "update Env...\n";
	for (auto spp : species_vec) spp->triggerPreCompute();
	env->computeEnv(t, this, S, dSdt);
}


// k = species_id
double Solver::calcSpeciesBirthFlux(int k, double t){
	if (debug) std::cout << "calc birthFlux...\n";
	auto spp = species_vec[k];	
	auto newborns_production = [this, spp](int i, double _t){
		double b1 = spp->birthRate(i, spp->getX(i), _t, env);
		return b1;	
	}; 
	double birthFlux = integrate_x(newborns_production, t, k);
	return birthFlux;	
}


vector<double> Solver::newborns_out(double t){  // TODO: make recompute env optional
	// update Environment from latest state
	//copyStateToCohorts(state.begin());  // not needed here because this is done in afterStep or upon cohorts update
	//env->computeEnv(t, this, state.begin(), rates.begin());
	updateEnv(t, state.begin(), rates.begin()); // this will trigger a precompute
//	for (int k = 0; k<species_vec.size(); ++k) preComputeSpecies(k,t);	

	vector<double> b_out;
	for (int k=0; k<species_vec.size(); ++k){	
		// calculate birthflux
		// [solved] wudx doesnt work here. Why?? - works now, no idea why it was not working earlier!
		//auto newborns_production = [this, k](int i, double t){
			//double z = species_vec[k]->getX(i);
			////species_vec[k]->preCompute(i,t,env);	
			//double b = species_vec[k]->birthRate(i,z,t,env);	
			////cout << "newborns of " << i << " = " << b << "\n"; 
			//return b; 
		//}; 
		//double birthFlux = integrate_x(newborns_production, t, k);
		//double birthFlux = integrate_wudx_above(newborns_production, t, 0, k);
		double birthFlux = calcSpeciesBirthFlux(k, t);
		b_out.push_back(birthFlux);
	}
	return b_out;
}

// FOR DEBUG ONLY, using TESTMODEL
vector<double> Solver::u0_out(double t){
	vector <double> u0out;
	vector <double> newbornsout = newborns_out(t);
	for (int k=0; k < species_vec.size(); ++k){
		u0out.push_back(newbornsout[k]/species_vec[k]->growthRate(-1, species_vec[k]->xb, t, env));
	}
	return u0out;
}



// xb
// |---*--|--*--|----*----|---*---|
//     0     1       2        3        <--- points
// 0      1     2         3       4    <--- breaks

// Test this function with this R script:
/**
x = c(0, 1,1.1,1.2, 2,2.2,2.3, 4.5,5.5,6.5)[length(x):1]
N = c(1, 1,1,  1,   2,2,  2,   3,   4,  5)[length(N):1]
plot(N~x)
abline(v=breaks, col="grey")
breaks = seq(0,10,1)
dens = rep(0, length(x))
J = length(x)

current_interval = length(breaks)-1
for (i in 1:J){
  while(breaks[current_interval] > x[i]) current_interval = current_interval-1
  dens[current_interval] = dens[current_interval]+N[i]
}
**/
struct point{
	double xmean = 0;
	double abund = 0;
	int    count = 0;
};	

std::vector<double> Solver::getDensitySpecies(int k, vector<double> breaks){
	auto spp = species_vec[k];
	vector <double> xx, uu;

	if (method == SOLVER_EBT){
		double xm = spp->getX(0)+1e-6;

		vector<point> points(breaks.size()-1);

		// assuming breaks are sorted ascending
		// cohorts are sorted descending
		int current_interval = breaks.size()-2;
		for (int i=0; i<spp->J; ++i){ // loop over all cohorts except boundary cohort
			double x = spp->getX(i);
			double N = spp->getU(i);
			if (i == spp->J-1) x = spp->xb + x/(N+1e-12); // real-ize x if it's boundary cohort
			
			while(breaks[current_interval]>x) --current_interval; // decrement interval until x fits
			//cout << current_interval << ", x = " << x << "(" << N << ") in [" << breaks[current_interval] << ", " << breaks[current_interval+1] << "]\n"; cout.flush();
			if (N>0){
				points[current_interval].abund += N;
				points[current_interval].count += 1;
				points[current_interval].xmean += N*x;
			}
		}

		
		for (int i=0; i<points.size(); ++i) if (points[i].count>0) points[i].xmean /= points[i].abund;
		
		// remove 0-count points
		auto pred = [this](const point& p) -> bool {
			return p.count == 0; 
		};

		auto pend = std::remove_if(points.begin(), points.end(), pred);
		points.erase(pend, points.end());

		//cout << "mean x and abund (removed):\n";
		//for (int i=0; i<points.size(); ++i) cout << i << "\t" << points[i].count << "\t" << points[i].xmean << "\t" << points[i].abund << "\n";	
		//cout << "--\n";

		if (points.size() > 2){
			vector<double> h(points.size());
			h[0] = (points[1].xmean+points[0].xmean)/2 - spp->xb;
			for (int i=1; i<h.size()-1; ++i) h[i] = (points[i+1].xmean - points[i-1].xmean)/2;
			h[h.size()-1] = xm - (points[h.size()-1].xmean+points[h.size()-2].xmean)/2;

			xx.reserve(points.size());
			uu.reserve(points.size());
			for (int i=0; i<points.size(); ++i){
				xx.push_back(points[i].xmean);
				uu.push_back(points[i].abund / h[i]);
			}
		}
	}

	else if (method == SOLVER_CM){
		xx.reserve(spp->J);
		uu.reserve(spp->J);
		for (int i=spp->J-1; i>=0; --i){
			xx.push_back(spp->getX(i));
			uu.push_back(spp->getU(i));
		}
		
	}	
	
	else {
		xx.reserve(spp->J);
		uu.reserve(spp->J);
		for (int i=0; i<spp->J-1; ++i){
			xx.push_back(spp->getX(i));
			uu.push_back(spp->getU(i));
		}
		
	}	
//	else {
//		throw std::runtime_error("This function can only be called for the EBT solver");
//	}

	if (xx.size() >= 2){ 
		
		Spline spl;
		spl.splineType = Spline::LINEAR; //Spline::CONSTRAINED_CUBIC;
		spl.extrapolate = Spline::QUADRATIC; //Spline::ZERO;
		spl.set_points(xx, uu);
		 
		vector <double> dens;
		dens.reserve(xx.size());
		for (int i=0; i<breaks.size(); ++i){
			dens.push_back(spl.eval(breaks[i]));			
		}
		
		return dens;
	}
	else {
		return vector<double>(breaks.size(), 0);
	}
}





