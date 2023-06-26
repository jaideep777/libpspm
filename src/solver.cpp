#include "solver.h"

#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include <algorithm>
#include <functional>
#include <fstream>
#include "util.h"

using namespace std;

// ~~~~~~~~~~~ SOLVER ~~~~~~~~~~~~~~~~~~~~~

std::map<std::string, PSPM_SolverType> Solver::methods_map = 
		{{"FMU",  SOLVER_FMU}, 
		 {"MMU",  SOLVER_MMU}, 
		 {"CM",   SOLVER_CM}, 
		 {"EBT",  SOLVER_EBT}, 
		 {"IFMU", SOLVER_IFMU}, 
		 {"ABM",  SOLVER_ABM}, 
		 {"IEBT", SOLVER_IEBT},
		 {"ICM",  SOLVER_ICM}};


Solver::Solver(PSPM_SolverType _method, string ode_method) : odeStepper(ode_method, 0, 1e-6, 1e-6) {
	method = _method;

	// FMU and ABM have only 1 internal state variable (x), rest have 2 (x,u)
	bool cond = (method == SOLVER_FMU || method == SOLVER_IFMU || method == SOLVER_ABM);
	n_statevars_internal = (cond)? 1:2;

}

Solver::Solver(std::string _method, std::string ode_method) : Solver(methods_map.at(_method), ode_method){
}


/// @brief Add the given species to the solver. 
/// @param xbreaks             breaks to use for discretization of the size axis
/// @param s                   species to add
/// @param n_extra_vars        number of cummulative variables to add for this species
/// @param input_birth_flux    The initial input birth flux for the species
/// @details This function creates metadata associated with the discretized size axis according to the specified solver, such as size at birth, initial number of cohorts/bins, and allocates space for the species in the state vector. 
void Solver::addSpecies(std::vector<double> xbreaks, Species_Base* s, int n_extra_vars, double input_birth_flux){
	s->set_inputBirthFlux(input_birth_flux);
	s->n_extra_statevars = n_extra_vars;

	if (method == SOLVER_FMU || method == SOLVER_IFMU || method == SOLVER_ABM){
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
	else if (method == SOLVER_ICM ) J = xbreaks.size();
	else if (method == SOLVER_EBT)  J = xbreaks.size();
	else if (method == SOLVER_IEBT) J = xbreaks.size();
	else if (method == SOLVER_ABM)  J = xbreaks.size();    // For ABM solver, this is a temporary size thats used to generate the initial density distribution. s will be resized during init to abm_n0.
	else    throw std::runtime_error("Unsupported method");

	s->x = xbreaks;
	s->resize(J);
	
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

	initializeSpecies(s);
}


void Solver::addSpecies(int _J, double _xb, double _xm, bool log_breaks, Species_Base* s,
							   int n_extra_vars, double input_birth_flux){
	vector<double> xbreaks;
	if (log_breaks) xbreaks = logseq(_xb, _xm, _J+1);
	else            xbreaks =    seq(_xb, _xm, _J+1);

	addSpecies(xbreaks, s, n_extra_vars, input_birth_flux);
}


void Solver::removeSpecies(Species_Base * spp){
	std::vector<Species_Base*>::iterator it = std::find(species_vec.begin(), species_vec.end(), spp);
	if (it != species_vec.end()){
		std::cout << "Removing species: " << spp << "\n";
		// Not freeing memory here: allocation of memory for species is done by user, so freeing should also be done by user
		species_vec.erase(it);
		resizeStateFromSpecies();
	}
	else{
		std::cout << "Species " << spp << " not found in the solver.\n";
	}
}


void Solver::addSystemVariables(int _s){
	n_statevars_system = _s;
	state.resize(state.size() + _s);
	rates.resize(state.size() + _s);
}


void Solver::resetState(double t0){  // FIXME: This is currently redundant, and needs to be improved with reset of both state and cohorts for a true reset of state
	current_time = t0;
	odeStepper.reset(t0, control.ode_eps, control.ode_eps); // = RKCK45<vector<double>> (0, control.ode_eps, control.ode_initial_step_size);  // this is a cheap operation, but this will empty the internal containers, which will then be (automatically) resized at next 1st ODE step. Maybe add a reset function to the ODE stepper? 

	// set birth time for each cohort to current_time
	for (auto s : species_vec){
		for (int i=0; i<s->J; ++i) s->set_birthTime(i, current_time); // FIXME: doesnt make sense, because larger cohorts would have been born earlier, but birthTime is not used anyways
	}

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
	if (method == SOLVER_ABM){
		for (auto spp : species_vec){
			for (int i=0; i<spp->J; ++i) maxsize = std::max(maxsize, spp->getX(i));
		}
	}
	else {
		for (auto spp : species_vec) maxsize = std::max(maxsize, spp->get_maxSize());
	}
	return maxsize;
}




void Solver::print(){
	std::cout << ">> SOLVER \n";
	string types[] = {"FMU", "MMU", "CM", "EBT", "Implicit FMU", "ABM", "Implicit EBT", "Implicit CM"};
	std::cout << "+ Type: " << types[method] << std::endl;

	std::cout << "+ State size = " << state.size() << "\n";
	std::cout << "+ Rates size = " << rates.size() << "\n";
	std::cout << "+ Environment = " << env << "\n";
	std::cout << "+ Species (" << species_vec.size() << "):\n";
	for (int i=0; i<species_vec.size(); ++i) {
		std::cout << "Sp (" << i << "):\n";
		species_vec[i]->print();
	}
	std::cout.flush();
}


/// @brief      initializes species - sets initial state (x, u, abc) for all cohorts depending on the solver
/// @param s    species to be initialized
// TODO: should this take t0 as an argument, instead of setting to 0? 
// [resolved] todo: maybe make use of copyCohortsToState here instead ofnmanually updating state elements?
void Solver::initializeSpecies(Species_Base * s){
		if (env == nullptr) throw std::runtime_error("Error: Environment has not been set. You must set it before addSpecies().");

		// set x and u of boundary cohort
		// Boundary cohort is not in state, but used as a reference.	
		s->set_xb(s->xb); // set x of boundary cohort 
		s->set_ub(0);     // set initial density of boundary cohort to 0.
		
		// set birth time for each cohort to current_time
		// FIXME: current_time has never been initialized till this point. It is only init in resetState() 
		for (int i=0; i<s->J; ++i) s->set_birthTime(i, current_time); // FIXME: doesnt make sense, because larger cohorts would have been born earlier, but birthTime is not used anyways

		// set x, u for all cohorts
		if (method == SOLVER_FMU || method == SOLVER_IFMU){
			for (size_t i=0; i<s->J; ++i){
				double X = s->X[i];			// for FMU, X is the midpoint of the cell edges
				double U = s->init_density(i, X, env); 
				s->setX(i,X); 
				s->setU(i,U);
				// *it++ = U;		// u in state (only)
			}
		}
		
		if (method == SOLVER_CM || method == SOLVER_ICM){
			for (size_t i=0; i<s->J; ++i){
				double X = s->x[i];
				double U = s->init_density(i, X, env); 
				s->setX(i,X); 
				s->setU(i,U);
				// *it++ = X;									// x in state
				// *it++ = (use_log_densities)? log(U) : U;	// u in state 
			}
		}
		
		if (method == SOLVER_EBT || method == SOLVER_IEBT){
			// x, u for internal cohorts 
			for (size_t i=0; i<s->J-1; ++i){
				double X = (s->x[i]+s->x[i+1])/2.0;			
				double U = s->init_density(i, X, env)*(s->x[i]-s->x[i+1]); 
				s->setX(i,X); 
				s->setU(i,U);
				// *it++ = X;	// x in state
				// *it++ = U;	// u in state 
			}
			// set pi0, N0 as x, u for the last cohort. This scheme allows using this last cohort with xb+pi0 in integrals etc 
			// *it++ = 0; *it++ = 0;
			s->setX(s->J-1,0); // [resolved] todo: should this be set to xb for init_state and set to 0 again later? maybe not, as init_state is not expected to be x dependent
			s->setU(s->J-1,0); 		
		}
		
		// FIXME: abm_n0 and x.size() can be different - we want higher x.size() to get a high-res initial density function, but when we draw from it, we only draw abm_n0 individuals
		if (method == SOLVER_ABM){
			// Create the initial density distribution from which we will draw individuals
			vector<double> Uvec;
			Uvec.reserve(s->x.size()-1);
			for (size_t i=0; i<s->x.size()-1; ++i){
				double X = s->x[i];
				cout << "i/X = " << i << "/" << X << endl;
				double U = s->init_density(i, X, env)*(s->x[i+1]-s->x[i]); 
				Uvec.push_back(U);	
			}
			//cout << "HERE\n";
			//for (size_t i=0; i<s->x.size()-1; ++i) cout << s->x[i] << " " << Uvec[i] << "\n";
		
			s->resize(control.abm_n0); // Once initial density dist has been obtained, resize species to n0

			std::discrete_distribution<int> distribution(Uvec.begin(), Uvec.end()); // for drawing intervals
			std::uniform_real_distribution<> distribution2(0,1);                         // for drawing X within interval

			// Utot = sum(Uvec) = sum(u[i] * dx[i])
			double Utot = std::accumulate(Uvec.begin(), Uvec.end(), 0.0, std::plus<double>());			
			double N_cohort = Utot/s->J;
			s->set_ub(N_cohort);
			for (int i=0; i<s->J; ++i){
				int id = distribution(generator);
				//cout << id << " ";
				double xi = s->x[id] + distribution2(generator)*(s->x[id+1]-s->x[id]);
				//cout << xi << "\n";
				s->setX(i, xi);
				s->setU(i, N_cohort);
			}

			s->sortCohortsDescending();
			
//			ofstream fout("abm_init.txt");
//			for (int i=0; i<s->J; ++i){
//				fout << s->getX(i) << "\t" << s->getU(i) << "\n";
//			}		
//			fout.close();	
		}

		// initialize extra state for each cohort ... previously also copied to state
		// For EBT, the initialization of extra variables may need state info, so need to realize pi0-cohort
		if (method == SOLVER_EBT || method == SOLVER_IEBT) realizeEbtBoundaryCohort(s);
		s->initExtraState(current_time, env);
		if (method == SOLVER_EBT || method == SOLVER_IEBT) restoreEbtBoundaryCohort(s);
		//if (s->n_extra_statevars > 0){  // FIXME: maybe redundant
			//auto it_prev = it;
			//s->init_ExtraState(it);  // this also inits the extra state of boundary cohort, but without advancing the iterator
			//assert(distance(it_prev, it) == s->n_extra_statevars*s->J); 
		//}

		resizeStateFromSpecies();
		copyCohortsToState();

}


void Solver::initialize(){
	// FIXME: Where is initialization of system vars?
//	vector<double>::iterator it = state.begin() + n_statevars_system; // TODO: replace with init_sState() 
	// for (int k=0; k<species_vec.size(); ++k){
	// 	Species_Base* s = species_vec[k];
	// 	initializeSpecies(s);
	// }
	resizeStateFromSpecies();
	copyCohortsToState();
}


void Solver::realizeEbtBoundaryCohort(Species_Base * spp){
	// backup pi0, N0 from last (youngest) cohort <-- cohorts are sorted descending
	pi0 = spp->getX(spp->J-1);
	N0  = spp->getU(spp->J-1);

	// real-ize pi0-cohort with actual x0 value
	double x0 = spp->xb + pi0/(N0+1e-12);
	spp->setX(spp->J-1, x0);
}


void Solver::restoreEbtBoundaryCohort(Species_Base * spp){
	// Copy saved value of pi0 back to the pi0-cohort (pi0 cohort is at index J-1)
	spp->setX(spp->J-1, pi0);
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
		if (method == SOLVER_CM || method == SOLVER_ICM){
			for (size_t i=0; i<s->J; ++i){
				double X = *it++;	// get x from state
				double U = *it++;	// get u from state
			    U = (use_log_densities)? exp(U) : U;	// u in cohorts 
				s->setX(i,X); 
				s->setU(i,U);
			}
		}
		if (method == SOLVER_EBT || method == SOLVER_IEBT){
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
		if (method == SOLVER_CM || method == SOLVER_ICM){
			for (size_t i=0; i<s->J; ++i){
				double X = s->getX(i); 
				double U = s->getU(i);
			    U = (use_log_densities)? log(U) : U;	// log(u) in state 
				*it++ = X;	// set x to state
				*it++ = U;	// set u to state
			}
		}
		if (method == SOLVER_EBT || method == SOLVER_IEBT){
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


void Solver::calcOdeRatesImplicit(double t, vector<double>::iterator S, vector<double>::iterator dSdt){
	vector<double>::iterator its = S    + n_statevars_system; // Skip system variables
	vector<double>::iterator itr = dSdt + n_statevars_system;
	
	for (int s = 0; s<species_vec.size(); ++s){
		auto spp = species_vec[s];
		
		its += (n_statevars_internal)*spp->J; // skip x and u 
		for (int i=0; i<n_statevars_internal*spp->J; ++i) *itr++ = 0; // set dx/dt and du/dt to 0 
	
		// FIXME: Maybe a good idea to realize pi0 cohort before calc extra rates
		if (spp->n_extra_statevars > 0){
			auto itr_prev = itr;
			spp->getExtraRates(itr); // TODO/FIXME: Does calc of extra rates need t and env?
			assert(distance(itr_prev, itr) == spp->n_extra_statevars*spp->J);
			its += spp->n_extra_statevars*spp->J; 	
		}
	}

}


void Solver::step_to(double tstop){
	auto func = [](double t){};
	step_to(tstop, func);
}


void Solver::updateEnv(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt){
	if (debug) std::cout << "update Env..." << std::endl;
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

std::vector<double> Solver::getDensitySpecies(int k, vector<double> breaks, Spline::Extr extrapolation_method){
	auto spp = species_vec[k];
	vector <double> xx, uu;
	
	if (method == SOLVER_ABM) spp->sortCohortsDescending(); // ABM needs sorted cohorts

	if (method == SOLVER_EBT || method == SOLVER_IEBT || method == SOLVER_ABM){ // EBT and ABM have very simular structure so use the same density calculation algo
		double xm = spp->getX(0)+1e-6;

		vector<point> points(breaks.size()-1);

		// assuming breaks are sorted ascending
		// and cohorts are sorted descending
		int current_interval = breaks.size()-2;
		for (int i=0; i<spp->J; ++i){ // loop over all cohorts except boundary cohort
			double x = spp->getX(i);
			double N = spp->getU(i);
			if (method == SOLVER_EBT) if (i == spp->J-1) x = spp->xb + x/(N+1e-12); // For EBT, real-ize x if it's boundary cohort
			
			while(breaks[current_interval]>x) --current_interval; // decrement interval until x fits
			//cout << current_interval << ", x = " << x << "(" << N << ") in [" << breaks[current_interval] << ", " << breaks[current_interval+1] << "]\n"; cout.flush();
			if (N>0){
				points[current_interval].abund += N;
				points[current_interval].count += 1;
				points[current_interval].xmean += N*x;
			}
		}

		// Compute mean x in each interval (each point corresponds to 1 interval)
		for (int i=0; i<points.size(); ++i) if (points[i].count>0) points[i].xmean /= points[i].abund;
		
		// remove 0-count points (i.e., delete intervals with no cohorts)
		auto pred = [this](const point& p) -> bool {
			return p.count == 0; 
		};

		auto pend = std::remove_if(points.begin(), points.end(), pred);
		points.erase(pend, points.end());

		//cout << "mean x and abund (removed):\n";
		//for (int i=0; i<points.size(); ++i) cout << i << "\t" << points[i].count << "\t" << points[i].xmean << "\t" << points[i].abund << "\n";	
		//cout << "--\n";

		// Now treat xmean as the x vector and calculate the width spanned by each point
		// to get u, divide abundance by width for each point
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

	else if (method == SOLVER_CM || method == SOLVER_ICM){
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

	// interpolate density at each value in breaks from uu
	if (xx.size() >= 2){ 
		
		Spline spl;
		spl.splineType = Spline::LINEAR; //Spline::CONSTRAINED_CUBIC;
		spl.extrapolate = extrapolation_method; //Spline::ZERO; //Spline::ZERO;
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


void Solver::save(std::ofstream &fout){
	fout << "Solver::v1\n";
	int m = static_cast<int>(method);
	fout << std::make_tuple(
		  m
		, n_statevars_internal
		, n_statevars_system
		, current_time
		, use_log_densities
		, pi0
		, N0);
	fout << '\n';

	fout << species_vec.size() << '\n';
	for (auto spp : species_vec) spp->save(fout);

	// we actually dont need the full state vector, as it can be reconstructed from cohorts. We only need the system variables
	fout << state; 

	odeStepper.save(fout);
}


void Solver::restore(std::ifstream &fin, vector<Species_Base*> spp_proto){
	string s; fin >> s; // version number (discard)
	assert(s == "Solver::v1");
	int m;
	fin >> m
	    >> n_statevars_internal
	    >> n_statevars_system
		>> current_time
		>> use_log_densities
		>> pi0
		>> N0;
	method = PSPM_SolverType(m);

	int nspecies;
	fin >> nspecies;
	assert(nspecies == spp_proto.size());
	species_vec = spp_proto;

	for (auto spp : species_vec){
		spp->restore(fin);
	}
	resizeStateFromSpecies(); // this includes system vars
	
	fin >> state;
	copyCohortsToState(); // this will overwrite species state (except system vars), but I guess this is better for sake of consistency
	
	odeStepper.restore(fin);
}

