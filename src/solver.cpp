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

//template<class Model, class Environment>
//const int Solver<Model,Environment>::xsize(){
//    if (method == SOLVER_FMU) return J;	
//    if (method == SOLVER_MMU) return J;  
//    if (method == SOLVER_CM ) return J;
//    if (method == SOLVER_EBT) return J;
//}

Solver::Solver(PSPM_SolverType _method, string ode_method) : odeStepper(ode_method, 0, 1e-6, 1e-6) {
	method = _method;
	//control.ode_method = ode_method;
}


////int Solver<Model,Environment>::setupLayout(Species<Model> &s){
//    // Set up layout of state vector, for e.g.
//    //  ------------------------------------------------------------
//    // | x | x | x : u | u | u : a | b | c | a | b | c | a | b | c |
//    //  ------------------------------------------------------------
//    //  In above layout, the internal variables x and u are tightly
//    //  packed, and the extra variables a, b, c are interleaved. 
//    //  This arrangement is cache friendly because:
//    //  1. When calculating integrals, x and u are traversed
//    //  2. When setting rates, typically a,b,c are calculated one 
//    //  after the other for each x.
//    s.clearVars();
//    if (method == SOLVER_FMU){ 
//        s.addVar("u", 1, s.J);		// uuuuu..
//    }
//    else{
//        s.addVar("X", 1, s.J);		// xxxxx.. uuuuu..
//        s.addVar("u", 1, s.J);		// 
//    }
//    for (int i=0; i < s.varnames_extra.size(); ++i){
//        s.addVar(s.varnames_extra[i], s.varnames_extra.size(), 1);  // abc abc abc .. 
//    }

//    return s.size();
//}


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

	// state.resize(state_size);   // state will be resized by addSpecies
	//rates.resize(state.size());
	
	
	std::fill(state.begin(), state.end(), 0); 
	std::fill(rates.begin(), rates.end(), -999); // DEBUG

	// initialize grid/cohorts for each species
	//for (int k=0; k<species_vec.size(); ++k){
		//Species_Base* s = species_vec[k]; // easy reference 

		//if (method == SOLVER_FMU){	
			//s->X.resize(s->J);
			//for (size_t i=0; i<s->J; ++i) s->X[i] = (s->x[i]+s->x[i+1])/2.0;
			
			//s->h.resize(s->J);	// This will be used only by FMU
			//for (size_t i=0; i<s->J; ++i) s->h[i] = s->x[i+1] - s->x[i];	
		//}

	//}
	//u0_out_history.clear();
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


//double Solver::get_u0(double t, int s){
	//Species_Base * spp = species_vec[s];
	
	//if (spp->bfin_is_u0in){
		//return spp->birth_flux_in;
	//}
	//else {	
		////double g = spp.mod->growthRate(spp.xb, t, env); // TODO: Move this computation to species
		////double u0 = (g>0)? spp.birth_flux_in * spp.mod->establishmentProbability(t, env)/g  :  0; 
		////return u0;
		//return 0;
	//}
//}


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
	//IteratorSet<vector<double>::iterator> iset(state.begin(), varnames.size(), vars, xsize());
	
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
	
		// set x for boundary cohort (BC is not in state, but used as a reference).	
		s->set_xb(s->xb);
		
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

	std::vector<double>::iterator it = state_begin + n_statevars_system; // no need to copy system state
	
	for (int k=0; k<species_vec.size(); ++k){
		Species_Base* s = species_vec[k];
		
		if (method == SOLVER_FMU || method == SOLVER_IFMU){
			for (size_t i=0; i<s->J; ++i){
				s->setU(i,*it++);
			}
		}
		if (method == SOLVER_CM){
			for (size_t i=0; i<s->J; ++i){
				double X = *it++;	// get x from state
				double U = *it++;	// get u from state
			    U = (use_log_densities)? exp(U) : U;	// u in state 
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
			    U = (use_log_densities)? log(U) : U;	// u in state 
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


////void Solver::getRatesFromCohorts(){

	//for (int k=0; k<species_vec.size(); ++k){
		//Species_Base* s = species_vec[k];
		
		//vector<double>::iterator it = rates.begin() + s->start_index;
		//if (method == SOLVER_FMU){
			//for (size_t i=0; i<s->J; ++i){
				//*it++ = s->getU(i);
			//}
		//}
		//if (method == SOLVER_CM){
			//for (size_t i=0; i<s->J; ++i){
				//double X = s->getX(i); 
				//double U = s->getU(i);
				//U = (use_log_densities)? exp(U) : U;	// u in state 
				//*it++ = X;
				//*it++ = U;
			//}
		//}
		//if (method == SOLVER_EBT){
			//// x, u for boundary and internal cohorts
			//for (size_t i=0; i<s->J; ++i){
				//double X = *it++; 
				//double U = *it++;
				//s->setX(i,X); 
				//s->setU(i,U);
			//}
		//}

		//if (s->n_extra_statevars > 0){  // If extra state variables have been requested, initialize them
			//auto it_prev = it;
			//s->copyExtraStateToCohorts(it);
			//assert(distance(it_prev, it) == s->n_extra_statevars*s->J); 
		//}
	//}
//}




////template<class Model, class Environment>
////void Solver<Model,Environment>::calcRates_extra(double t, vector<double>&S, vector<double>& dSdt){
////    //auto is = createIterators_state(S);
////    //auto ir = createIterators_rates(dSdt);
////    //auto& itx = is.get("X");
////    //auto& itre = ir.get(varnames_extra[0]);
	
////    //for (is.begin(), ir.begin(); !is.end(); ++is, ++ir){
////    //    auto it_returned = mod->calcRates_extra(t, *itx, itre);
////    //    assert(distance(itre, it_returned) == varnames_extra.size());
////    //}
////}


void Solver::step_to(double tstop){
	auto func = [](double t){};
	step_to(tstop, func);
}



void Solver::preComputeSpecies(int k, double t){
	auto spp = species_vec[k];

	if (method == SOLVER_EBT){ // for EBT, we need to pi0-cohort too.
		// backup and real-ize pi0-cohort
		// get pi0, N0 from last cohort
		double pi0  =  spp->getX(spp->J-1);
		double N0   =  spp->getU(spp->J-1);
		
		// update pi0-cohort with actual x0 value
		double x0 = spp->xb + pi0/(N0+1e-12);
		spp->setX(spp->J-1, x0);
	
		// precompute cohorts
		spp->preComputeAllCohorts(t,env);

		// restore pi0-cohort
		spp->setX(spp->J-1, pi0);
	}
	else{
		spp->preComputeAllCohorts(t,env);
	}
	
}


// k = species_id
double Solver::calcSpeciesBirthFlux(int k, double t){
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
	//copyStateToCohorts(state.begin());
	env->computeEnv(t, this, state.begin(), rates.begin());
	for (int k = 0; k<species_vec.size(); ++k) preComputeSpecies(k,t);	

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
		//species_vec[k]->preCompute(-1, t, env); // not req because precomputeAllCohorts called in newborns_out() precomputes BC too.	
		u0out.push_back(newbornsout[k]/species_vec[k]->growthRate(-1, species_vec[k]->xb, t, env));
	}
	return u0out;
}

//[>
//template<class Model, class Environment>
//double Solver<Model,Environment>::get_u0_out(){
//    return u0_out_history.back();
//}


//template<class Model, class Environment>
//double Solver<Model,Environment>::stepToEquilibrium(){
//    for (double t=0.05; ; t=t+0.05) {
//        step_to(t);
		
//        u0_out_history.push_back(u0_out());
//        if (u0_out_history.size() > 5) u0_out_history.pop_front();
//        //cout << t << " | ";  
//        //for (auto u : u0_out_history) cout << u << " "; // cout << "\n";
//        double max_err = -1e20;
//        //cout << u0_out_history.size() << endl; 
//        for (auto it = ++u0_out_history.begin(); it != u0_out_history.end(); ++it){
//            double err = *it - *(prev(it)); 
//            max_err = std::max(max_err, abs(err));
//            //cout << err << " ";
//        }
//        //cout << " | " << max_err << endl;
		
//        if (abs(max_err) < control.convergence_eps){
//            cout << t << " | " << J << " | "; for (auto u : u0_out_history) cout << u << " "; cout << "|" << max_err << endl;
//            break;
//        }	
//    }
	
//}



//template<class Model, class Environment>
//double Solver<Model,Environment>::createSizeStructuredVariables(vector<std::string> names){
//    varnames_extra = names;
//    resetState(x);
//}


//*/



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

std::vector<double> Solver::getDensitySpecies_EBT(int k, vector<double> breaks){
	auto spp = species_vec[k];

	//cout << "HRER" << endl;
	
	if (method == SOLVER_EBT){
		double xm = spp->getX(0)+1e-6;

		//vector<double> breaks = (logscale)?  logseq(spp->xb, xm, nbreaks) : seq(spp->xb, xm, nbreaks);
		
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
			return p.count == 0; // TODO: make conatiner-type safe
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

			vector <double> xx, uu;
			xx.reserve(points.size());
			uu.reserve(points.size());
			for (int i=0; i<points.size(); ++i){
				xx.push_back(points[i].xmean);
				uu.push_back(points[i].abund / h[i]);
			}
			
			Spline spl;
			spl.splineType = Spline::LINEAR; //Spline::CONSTRAINED_CUBIC;
			spl.extrapolate = Spline::ZERO;
			spl.set_points(xx, uu);
			
			vector <double> dens;
			dens.reserve(points.size());
			for (int i=0; i<breaks.size(); ++i){
				dens.push_back(spl.eval(breaks[i]));			
			}
			
			return dens;
		}
		else return vector<double>(breaks.size(), 0);
	}
	else {
		throw std::runtime_error("This function can only be called for the EBT solver");
	}

}





