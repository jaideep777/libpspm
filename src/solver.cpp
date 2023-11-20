#include "solver.h"
#include "index_utils.h"
#include "mcmc.h"

#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include <algorithm>
#include <functional>
#include <fstream>

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

inline std::vector <double> mids(const vector <double>& breaks){
	std::vector<double> a(breaks.size()-1);
	for (size_t i=0; i<a.size(); ++i) a[i] = (breaks[i]+breaks[i+1])/2;
	return a;
}

inline std::vector <double> left_edge(const vector <double>& breaks){
	std::vector<double> a(breaks.size()-1);
	for (size_t i=0; i<a.size(); ++i) a[i] = breaks[i];
	return a;
}

inline std::vector <double> diff(const vector <double>& breaks){
	std::vector<double> a(breaks.size()-1);
	for (size_t i=0; i<a.size(); ++i) a[i] = breaks[i+1]-breaks[i];
	return a;
}



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

	// JJ: Commenting out because the variable seems useless now, will always be 1 (u), as x is determined later per species
	// // FMU and ABM have only 1 internal state variable (x), rest have 2 (x,u)
	// bool cond = (method == SOLVER_FMU || method == SOLVER_IFMU || method == SOLVER_ABM);
	// n_statevars_internal = (cond)? 0:1; // Now we'll consider this as extra state variables - FIXME JJ: This doesnt work, because you're setting this to 0 for FMU and 1 for EBT...  

}

Solver::Solver(std::string _method, std::string ode_method) : Solver(methods_map.at(_method), ode_method){
}


/// @brief Add the given species to the solver. 
/// @param xbreaks             breaks to use for discretization of each state axis, sorted ascending
/// @param s                   species to add
/// @param _n_accumulators     number of cummulative variables to add for this species
/// @param input_birth_flux    The initial input birth flux for the species
/// @details This function creates metadata associated with the discretized size axis according to the specified solver, such as size at birth, initial number of cohorts/bins, and allocates space for the species in the state vector. 
void Solver::addSpecies(std::vector<std::vector<double>> xbreaks, Species_Base* s, int _n_accumulators, double input_birth_flux){
	// Ensure that as many coordinates as state variables have been supplied
	assert(xbreaks.size() == s->istate_size);

	// JJ Note: nD CM must be disallowed, because its integral (as of now) needs ordering of cohorts.
	if (method == SOLVER_CM || method == SOLVER_ICM){
		if (s->istate_size > 1) throw std::runtime_error("With the CM and ICM solvers, only 1 state variable per species is allowed.");
	}

	// JJ Note: nD FMU must be disallowed, because its boundary condition still needs to be sorted out.
	if (method == SOLVER_FMU || method == SOLVER_IFMU){
		if (s->istate_size > 1) throw std::runtime_error("With grid-based solvers, currently only 1D state variables are supported.");
	}

	// TODO: Should we check that all coordinates are sorted ascending?
	s->clear_vectors();

	s->set_inputBirthFlux(input_birth_flux);
	s->n_accumulators = _n_accumulators;

	// Create grid centres and grid dx along each axis
	s->x = xbreaks;
	for (int i=0; i<xbreaks.size(); ++i){
		s->X.push_back(mids(xbreaks[i]));
		s->h.push_back(diff(xbreaks[i]));
	}

	// in IFMU solver, if gridcells are labelled by lower edge, set X by lower edge rather than centre
	if (method == SOLVER_IFMU && !control.ifmu_centered_grids){
		s->X.clear(); // clear stuff set above
		for (int i=0; i<xbreaks.size(); ++i){
			s->X.push_back(left_edge(xbreaks[i]));
		}
	} 


	// Set species birth size - xb
	s->xb.clear();
	if (method == SOLVER_IFMU && control.ifmu_centered_grids){ 
		// For IFMU with centered grids, xb is the centre of the first nD cell 
		for (int i=0; i<s->istate_size; ++i) s->xb.push_back(s->X[i][0]);
	}
	else{  
		// For all other solvers, xb is the lower edge of the first nD cell 
		for (int i=0; i<s->istate_size; ++i) s->xb.push_back(s->x[i][0]);
	}


	// Create cohorts
	// vector<int> dim_centres, dim_edges;
	for (auto& vx : xbreaks) s->dim_centres.push_back(vx.size()-1);  // For each dimension, N breaks give rise to N-1 cohorts/cells 
	for (auto& vx : xbreaks) s->dim_edges.push_back(vx.size());  // For each dimension, N breaks give rise to N cell edges 

	s->n_grid_centres = std::accumulate(s->dim_centres.begin(), s->dim_centres.end(), 1, std::multiplies<int>());
	s->n_grid_edges   = std::accumulate(s->dim_edges.begin(),   s->dim_edges.end(),   1, std::multiplies<int>());
	if (s->n_grid_edges > 1e6) cout << "**** WARNING ****: The number of cohorts/cells may exceed 1M. Consider using a lower resolution\n\n";

	std::cout << "Find J" << std::endl;
	int J = 1;
	if      (method == SOLVER_FMU)   J = s->n_grid_centres; // xbreaks.size()-1;	
	else if (method == SOLVER_IFMU)  J = s->n_grid_centres;	
	else if (method == SOLVER_MMU)   J = s->n_grid_centres;  
	else if (method == SOLVER_CM )   J = s->n_grid_edges;  // For CM, its not strictly necessary to use grid edges, but this helps with tests on the Plant Model 
	else if (method == SOLVER_ICM )  J = s->n_grid_edges;
	else if (method == SOLVER_EBT)   J = s->n_grid_centres+1;  // As many cohorts as grid centers + 1 boundary cohort
	else if (method == SOLVER_IEBT)  J = s->n_grid_centres+1;
	else if (method == SOLVER_ABM)   J = s->n_grid_centres;    // For ABM solver, this is a temporary size thats used to generate the initial density distribution. s will be resized during init to abm_n0. FIXME JJ: Can ABM init be kept identical to EBT?
	else    throw std::runtime_error("Unsupported method");

	std::cout << "Resize with J" << std::endl;
	s->resize(J);

	std::cout << "Add species to vector" << std::endl;
	species_vec.push_back(s);

	std::cout << "Initialise species" << std::endl;
	initializeSpecies(s);

	// Test Print out X, x and h from the new species

	// std::cout << "Test: Solver::addSpecies: Print X/x/h" << std::endl;

	// std::cout << "Test: Solver::addSpecies: Print X/x/h" << std::endl;
	// for (int i=0; i<s->X.size(); ++i) std::cout << s->X[i] << std::endl;
	// for (int i=0; i<s->X.size(); ++i) std::cout << s->x << std::endl;
	// for (int i=0; i<s->X.size(); ++i) std::cout << s->h << std::endl;

}


void Solver::addSpecies(std::vector<int> _J, std::vector<double> _xb, std::vector<double> _xm, std::vector<bool> log_breaks, Species_Base* s, 
								int _n_accumulators, double input_birth_flux){

	if (_xb.size() != s->istate_size){
		throw std::runtime_error("Error: \nSolver::addSpecies: number of elements in xb doesnt match the species istate_size"); // Fix this to be more informative
	}
	if (_xb.size() != _xm.size()){
		throw std::runtime_error("Error: \nSolver::addSpecies: size of lower boundary and upper boundary for states doesn't match"); // Fix this to be more informative
	}
	if (_xb.size() != log_breaks.size()){
		throw std::runtime_error("Error: \nSolver::addSpecies: size of lower boundary and break information log_breaks for states doesn't match"); // Fix this to be more informative
	}
	// if (_xb.size()>1 && this->method != SOLVER_EBTN){
	// 	throw std::runtime_error("Error: \nSolver::addSpecies: nD solver for method not supported"); // Fix this to be more informative
	// }
	// Initialise as a grid and assume each dimension starts off with _J numbers

	std::cout << "Initialise as a grid and assume each dimension starts off with _J[k] numbers" << std::endl;
	int total_states = std::accumulate(_J.begin(), _J.end(), 1, std::multiplies<int>()) +1;

	//  Initialise as a grid
	std::cout << "Initialise as a grid" << std::endl;
	std::vector<std::vector<double>> breaks;

	for (int k=0; k< s->istate_size; ++k){
		std::vector<double> ax_breaks(_J[k]+1);
		if (log_breaks[k]) {
			ax_breaks = logseq(_xb[k], _xm[k], _J[k] + 1);
		}
		else {
			ax_breaks = seq(_xb[k], _xm[k], _J[k] + 1);
		}
		breaks.push_back(ax_breaks);
	}


	// std::cout << "Populate Grid" << std::endl;
	// xnbreaks[0] = _xb;
	// for (int i=1; i < xnbreaks.size(); i++){
	// 	std::vector<int> ind = index(i-1, dims);
	// 	std::vector<double> xn;
	// 	for(int k = 0; k < ind.size(); k++){
	// 		xn.push_back(breaks[k][ind[k]+1]);
	// 	}
	// 	xnbreaks[i] = xn;
	// }

	std::cout << "Grid populated:" << std::endl;
	for (auto& b : breaks) std::cout << b << '\n';
	std::cout.flush();

	addSpecies(breaks, s, _n_accumulators, input_birth_flux);
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


void Solver::addSystemVariables(const std::vector<double>& s_state_0){
	s_state = s_state_0;
	n_statevars_system = s_state_0.size();
	resizeStateFromSpecies();
	copyCohortsToState();
}


// void Solver::resetState(double t0){  // FIXME: This is currently redundant, and needs to be improved with reset of both state and cohorts for a true reset of state
// 	current_time = t0;
// 	odeStepper.reset(t0, control.ode_eps, control.ode_eps); // = RKCK45<vector<double>> (0, control.ode_eps, control.ode_initial_step_size);  // this is a cheap operation, but this will empty the internal containers, which will then be (automatically) resized at next 1st ODE step. Maybe add a reset function to the ODE stepper? 

// 	// set birth time for each cohort to current_time
// 	for (auto s : species_vec){
// 		for (int i=0; i<s->J; ++i) s->set_birthTime(i, current_time); // FIXME: doesnt make sense, because larger cohorts would have been born earlier, but birthTime is not used anyways
// 	}

// 	std::fill(state.begin(), state.end(), 0); 
// 	std::fill(rates.begin(), rates.end(), -999); // DEBUG
// }


void Solver::resizeStateFromSpecies(){
	int state_size_new = n_statevars_system;
	for (auto& spp : species_vec){
		state_size_new += n_statevars(spp);
		state_size_new += spp->J*(spp->n_accumulators); // add accumulators for all solvers
	}
	
	// std::cout << "In Solver::resizeStateFromSpecies before resize " <<std::endl;
	state.resize(state_size_new, -999);
	rates.resize(state_size_new, -999);
	// std::cout << "In Solver::resizeStateFromSpecies after resize " <<std::endl;
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

int Solver::n_statevars(Species_Base * spp){
	int n = spp->J*(1); // add u for all solvers
	if (method != SOLVER_FMU && method != SOLVER_IFMU){
		n += spp->J*(spp->istate_size);    // add x for solvers other than FMU family
	}
	return n;
}

int Solver::n_statevars_cohort(Species_Base * spp){
	int n = 1; // add u for all solvers
	if (method != SOLVER_FMU && method != SOLVER_IFMU){
		n += spp->istate_size;    // add x for solvers other than FMU family
	}
	return n;
}


// double Solver::maxSize(){
// 	double maxsize = 0;
// 	if (method == SOLVER_ABM){
// 		for (auto spp : species_vec){
// 			for (int i=0; i<spp->J; ++i) maxsize = std::max(maxsize, spp->getX(i));
// 		}
// 	}
// 	else {
// 		for (auto spp : species_vec) maxsize = std::max(maxsize, spp->get_maxSize());
// 	}
// 	return maxsize;
// }




void Solver::print(){
	std::cout << ">> SOLVER \n";
	string types[] = {"FMU", "MMU", "CM", "EBT", "IFMU", "ABM", "IEBT", "ICM", "EBTN", "IEBTN"};
	std::cout << "+ Type: " << types[method] << std::endl;

	std::cout << "+ State size = " << state.size() << "\n";
	std::cout << "+ Rates size = " << rates.size() << "\n";
	std::cout << "+ Environment = " << env << "\n";
	std::cout << "+ S-State variables = " << s_state << "\n";
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
// [resolved] todo: maybe make use of copyCohortsToState here instead of manually updating state elements?
void Solver::initializeSpecies(Species_Base * s){
		if (env == nullptr) throw std::runtime_error("Error: Environment has not been set. You must set it before addSpecies().");

		// set x and u of boundary cohort
		// Boundary cohort is not in state, but used as a reference.	

		std::cout << "Set up boundary" << std::endl;
		s->set_xb(s->xb); // set x of boundary cohort 
		std::cout << "Set up u at boundary" << std::endl;
		s->set_ub(0);     // set initial density of boundary cohort to 0.
		
		std::cout << "set birthtime" << std::endl;
		// set birth time for each cohort to current_time
		// FIXME: current_time has never been initialized till this point. It is only init in resetState() 
		for (int i=0; i<s->J; ++i) s->set_birthTime(i, current_time); // FIXME: doesnt make sense, because larger cohorts would have been born earlier, but birthTime is not used anyways


		std::cout << "set x and u for all cohorts" << std::endl;
		// set x, u for all cohorts
		if (method == SOLVER_FMU || method == SOLVER_IFMU){
			for (size_t i=0; i<s->J; ++i){
				vector<double> X = id_utils::coord_value(id_utils::index(i, s->dim_centres), s->X);
				s->setX(i,X); 

				double U = s->init_density(i, env); 
				s->setU(i,U);
				// *it++ = U;		// u in state (only)
			}
		}
		
// 		if (method == SOLVER_CM || method == SOLVER_ICM){
// 			for (size_t i=0; i<s->J; ++i){
// 				double X = s->x[i];
// 				s->setX(i,X); 

// 				double U = s->init_density(i, X, env); 
// 				s->setU(i,U);
// 				// *it++ = X;									// x in state
// 				// *it++ = (use_log_densities)? log(U) : U;	// u in state 
// 			}
// 		}
		
		if (method == SOLVER_EBT || method == SOLVER_IEBT){
			// x, u for internal cohorts 
			for (size_t i=0; i<s->J-1; ++i){
				vector<double> X = id_utils::coord_value(id_utils::index(i, s->dim_centres), s->X);
				s->setX(i,X); 

				vector<double> dx = id_utils::coord_value(id_utils::index(i, s->dim_centres), s->h);
				cout << "dx = " << dx << '\n';
				double dV = std::accumulate(dx.begin(), dx.end(), 1.0, std::multiplies<double>());
				cout << "dV = " << dV << '\n';
				double U = s->init_density(i, env)*dV; 
				s->setU(i,U);
				cout << "Init: X = " << X << " / U = " << U << '\n';
			}
			// set pi0, N0 as x, u for the last cohort. This scheme allows using this last cohort with xb+pi0 in integrals etc 
			s->setX(s->J-1, vector<double>(s->istate_size, 0)); 
			s->setU(s->J-1, 0); 		
		}

		// FIXME: abm_n0 and x.size() can be different - we want higher x.size() to get a high-res initial density function, but when we draw from it, we only draw abm_n0 individuals
		// FIXME: Note that X must be set before calculating U.
		if (method == SOLVER_ABM){
			// Simulate initial elements using a Monte Carlo simulation 

			/* Find total density - cheat by using the initial grid a little  */
			double Utot = 0.0; 
			for (size_t i=0; i<s->J-1; ++i){
				vector<double> dx = id_utils::coord_value(id_utils::index(i, s->dim_centres), s->h);
				double dV = std::accumulate(dx.begin(), dx.end(), 1.0, std::multiplies<double>());
				double U = s->init_density(i, env)*dV; 
				Utot += U;
			}

			// Since we don't have an 'independent' density function and we will in any case substitute 
			// the x values. there should hopefully be a better way of doing this but since it is a one time
			// simulation it might be ok to keep it this way
			auto targetDensity = [s, this](const std::vector<double>& x) {
				s->setX(s->J-2, x);
				return s->init_density(s->J-2, env);
   	 		};

			// Create sampler chain

			std::vector<double> x_min = s->getX(s->J-1);
			std::vector<double> x_max = s->getX(0);
			std::vector<double> sd(x_min.size());
			for(size_t k = 0; k < x_min.size(); ++k){
				sd[k] = (x_max[k] - x_min[k])/10;
			}

			MCMCSampler sampler(x_min, x_max, sd, control.abm_numChains, control.abm_burnin, 1);
			sampler.run_chains(targetDensity, control.abm_n0 + control.abm_burnin);

			/* here instert test for convergence */

			/* allocate new samples */
			std::vector<std::vector<double>> sample_x = sampler.sample(control.abm_n0);
			s->resize(control.abm_n0 + 1); // Elisa: adding one for the boundary cohort - is this necessary?

			// allocate values to all elements
			for(size_t i = 0; i < sample_x.size(); ++i){
				s->setX(i, sample_x[i]);
			}
			//Create the initial density distribution from which we will draw individuals
			double N_cohort = Utot/s->J; 
			s->set_ub(N_cohort);
			for (int i=0; i<s->J; ++i){ 
				s->setU(i, N_cohort);
			}
			
//			ofstream fout("abm_init.txt");
//			for (int i=0; i<s->J; ++i){
//				fout << s->getX(i) << "\t" << s->getU(i) << "\n";
//			}		
//			fout.close();	
		}

		// initialize extra state for each cohort ... previously also copied to state
		// For EBT, the initialization of extra variables may need state info, so need to realize pi0-cohort
		if (method == SOLVER_EBT || method == SOLVER_IEBT) realizeEbtBoundaryCohort(s);
		s->initAccumulators(current_time, env);
		if (method == SOLVER_EBT || method == SOLVER_IEBT) restoreEbtBoundaryCohort(s);
		//if (s->n_accumulators > 0){  // FIXME: maybe redundant
			//auto it_prev = it;
			//s->init_ExtraState(it);  // this also inits the extra state of boundary cohort, but without advancing the iterator
			//assert(distance(it_prev, it) == s->n_accumulators*s->J); 
		//}

		resizeStateFromSpecies();
		copyCohortsToState();

}


// void Solver::initialize(){
// 	// [resolved FIXME]: Where is initialization of system vars? - Now addSystemVariables() takes initial values as argument
// //	vector<double>::iterator it = state.begin() + n_statevars_system; // TODO: replace with init_sState() 
// 	// for (int k=0; k<species_vec.size(); ++k){
// 	// 	Species_Base* s = species_vec[k];
// 	// 	initializeSpecies(s);
// 	// }
// 	resizeStateFromSpecies();
// 	copyCohortsToState();
// }


// //everything is just doubled for xn here
// void Solver::realizeEbtBoundaryCohort(Species_Base * spp){
// 	// backup pi0, N0 from last (youngest) cohort <-- cohorts are sorted descending
// 	pi0 = spp->getX(spp->J-1);
// 	N0  = spp->getU(spp->J-1);

// 	// real-ize pi0-cohort with actual x0 value
// 	double x0 = spp->xb + pi0/(N0+1e-12);
// 	spp->setX(spp->J-1, x0);
// }


void Solver::realizeEbtBoundaryCohort(Species_Base * spp){
	// backup pi0, N0 from last (youngest) cohort <-- cohorts are sorted descending
	pi0 = spp->getX(spp->J-1);
	N0  = spp->getU(spp->J-1);

	// real-ize pi0-cohort with actual x0 value
	std::vector<double> x0 = spp->xb;
	for(size_t k = 0; k < x0.size(); ++k){
		x0[k] += pi0[k]/(N0+1e-12);
	}
	spp->setX(spp->J-1, x0);
}


// void Solver::restoreEbtBoundaryCohort(Species_Base * spp){
// 	// Copy saved value of pi0 back to the pi0-cohort (pi0 cohort is at index J-1)
// 	spp->setX(spp->J-1, pi0);
// }

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

	std::vector<double>::iterator it = state_begin;
	
	// copy system state
	for (int i=0; i<n_statevars_system; ++i) s_state[i] = *it++; 

	size_t current_state = n_statevars_system;
	size_t state_size = state.size();
	// std::cout << "state size ---> " << state.size() << std::endl;
	// std::cout << "After state: " << n_statevars_system << std::endl;

	for (int k=0; k<species_vec.size(); ++k){
		Species_Base* s = species_vec[k];
		
		// std::cout << "In species " << k << std::endl;

		s->set_xb(s->xb); // Important: "touch" boundary cohort to trigger a precompute. Note that setSize() triggers it.
		if (method == SOLVER_FMU || method == SOLVER_IFMU){
			for (size_t i=0; i<s->J; ++i){
				s->setU(i,*it++);
			}
		}
		if (method == SOLVER_CM || method == SOLVER_ICM){
			for (size_t i=0; i<s->J; ++i){
				vector<double> X(s->istate_size);
				for (double& xx : X) xx = *it++;	// get x from state
				double U = *it++;	// get u from state
			    U = (use_log_densities)? exp(U) : U;	// u in cohorts 
				s->setX(i,X); 
				s->setU(i,U);
			}
		}
		if (method == SOLVER_EBT || method == SOLVER_IEBT){
			// x, u for boundary and internal cohorts
			for (size_t i=0; i<s->J; ++i){
				vector<double> X(s->istate_size);
				for (double& xx : X) xx = *it++;	// get x from state
				double U = *it++;
				s->setX(i,X); 
				s->setU(i,U);
			}
		}

		if (s->n_accumulators > 0){  // If extra state variables have been requested, initialize them
			auto it_prev = it;
			s->copyAccumulatorsToCohorts(it);
			assert(distance(it_prev, it) == s->n_accumulators*s->J); 
		}
	}
}


void Solver::copyCohortsToState(){
	if (debug) std::cout << "state <--- cohorts\n";
	// std::cout << "state size: " << state.size() << std::endl;
	std::vector<double>::iterator it = state.begin();
	
	// copy system state
	for (int i=0; i<n_statevars_system; ++i) *it++ = s_state[i]; 
	
	for (int k=0; k<species_vec.size(); ++k){
		Species_Base* s = species_vec[k];

		// std::cout << "Species size: " << s->J << std::endl;
		
		if (method == SOLVER_FMU || method == SOLVER_IFMU){
			for (size_t i=0; i<s->J; ++i){
				*it++ = s->getU(i);
			}
		}
		if (method == SOLVER_CM || method == SOLVER_ICM){
			for (size_t i=0; i<s->J; ++i){
				vector<double> X = s->getX(i); 
				double U = s->getU(i);
			    U = (use_log_densities)? log(U) : U;	// log(u) in state 
				for (double& xx : X) *it++ = xx;	// get x from state
				*it++ = U;	// set u to state
			}
		}
		if (method == SOLVER_EBT || method == SOLVER_IEBT){
			// x, u for boundary and internal cohorts
			for (size_t i=0; i<s->J; ++i){
				vector<double> X = s->getX(i); 
				double U = s->getU(i);
				for (double& xx : X) *it++ = xx;	// get x from state
				*it++ = U;
			}
		}

		if (s->n_accumulators > 0){  // If extra state variables have been requested, initialize them
			auto it_prev = it;
			s->copyAccumulatorsToState(it);
			assert(distance(it_prev, it) == s->n_accumulators*s->J); 
		}
	}
}


void Solver::calcOdeRatesImplicit(double t, vector<double>::iterator S, vector<double>::iterator dSdt){
	vector<double>::iterator its = S    + n_statevars_system; // Skip system variables
	vector<double>::iterator itr = dSdt + n_statevars_system;
	
	for (int s = 0; s<species_vec.size(); ++s){
		auto spp = species_vec[s];
		
		its += n_statevars(spp); // skip x and u 
		for (int i=0; i<n_statevars(spp); ++i) *itr++ = 0; // set dx/dt and du/dt to 0 
	
		// [resolved]: Maybe a good idea to realize pi0 cohort before calc extra rates
		if (spp->n_accumulators > 0){
			if (method == SOLVER_EBT || method == SOLVER_IEBT) realizeEbtBoundaryCohort(spp);
			auto itr_prev = itr;
			spp->accumulatorRates(itr); // TODO/FIXME: Does calc of extra rates need t and env?
			assert(distance(itr_prev, itr) == spp->n_accumulators*spp->J);
			its += spp->n_accumulators*spp->J; 	
			if (method == SOLVER_EBT || method == SOLVER_IEBT) restoreEbtBoundaryCohort(spp);	
		}
	}

}


void Solver::step_to(double tstop){
	auto func = [](double t){};
	step_to(tstop, func);
}


void Solver::updateEnv(double t, std::vector<double>::iterator S, std::vector<double>::iterator dSdt){
	//if (debug) std::cout << "update Env..." << std::endl;
	for (auto spp : species_vec) spp->triggerPreCompute();
	env->computeEnv(t, this, S, dSdt);
}


// k = species_id
double Solver::calcSpeciesBirthFlux(int k, double t){
	// if (debug) std::cout << "calc birthFlux...\n";
	auto spp = species_vec[k];	
	auto newborns_production = [this, spp](int i, double _t){
		double b1 = spp->birthRate(i, _t, env);
		return b1;	
	}; 
	double birthFlux = state_integral(newborns_production, t, k);
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
		u0out.push_back(newbornsout[k]/species_vec[k]->growthRate(-1, t, env)[0]);
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

// breaks must be sorted ascending
std::vector<double> Solver::getDensitySpecies1D(int k, int dim, const std::vector<double>& breaks, Spline::Extr extrapolation_method){
	auto spp = species_vec[k];
	vector <double> xx, uu;
	
	// if (method == SOLVER_ABM) spp->sortCohortsDescending(dim); // ABM needs sorted cohorts
	// if (method == SOLVER_EBT || method == SOLVER_IEBT) spp->sortCohortsDescending(dim, 1); // for EBT, skip boundary cohort and sort

	if (method == SOLVER_EBT || method == SOLVER_IEBT || method == SOLVER_ABM){ // EBT and ABM have very simular structure so use the same density calculation algo
		double xm = spp->getX(0)[dim]+1e-6;

		vector<point> points(breaks.size());

		// int current_interval = breaks.size()-2;
		for (int i=0; i<spp->J; ++i){ // loop over all cohorts 
			double x = spp->getX(i)[dim];
			double N = spp->getU(i);
			if (method == SOLVER_EBT || method == SOLVER_IEBT) if (i == spp->J-1) x = spp->xb[dim] + x/(N+1e-12); // For EBT, real-ize x if it's boundary cohort
			
			int current_interval = std::upper_bound(breaks.begin(), breaks.end(), x) - breaks.begin() -1; // this assumes breaks is sorted ascending
			current_interval = std::clamp<int>(current_interval, 0, breaks.size()-1);
			// while(breaks[current_interval]>x) --current_interval; // decrement interval until x fits
			// cout << current_interval << ", x = " << x << "(" << N << ") in [" << breaks[current_interval] << ", " << breaks[current_interval+1] << "]\n"; cout.flush();
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

		// cout << "mean x and abund (removed):\n";
		// for (int i=0; i<points.size(); ++i) cout << i << "\t" << points[i].count << "\t" << points[i].xmean << "\t" << points[i].abund << "\n";	
		// cout << "--\n";

		// Now treat xmean as the x vector and calculate the width spanned by each point, i.e., 1D Voronoi cells
		// to get u, divide abundance by width for each point
		if (points.size() > 2){
			vector<double> h(points.size());
			h[0] = (points[1].xmean+points[0].xmean)/2 - spp->xb[dim];
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
			xx.push_back(spp->getX(i)[dim]);
			uu.push_back(spp->getU(i));
		}
		
	}	
	
	else {
		xx.reserve(spp->J);
		uu.reserve(spp->J);
		for (int i=0; i<spp->J-1; ++i){
			xx.push_back(spp->getX(i)[dim]);
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


void Solver::save(std::ostream &fout){
	fout << "Solver::v2\n";
	int m = static_cast<int>(method);
	fout << std::make_tuple(
		  m
		, n_statevars_system
		, current_time
		, use_log_densities
		, N0);
	fout << '\n';
	fout << pi0 << '\n';
	// [resolved] we dont need the full state vector, as it can be reconstructed from cohorts. We only need the system variables
	fout << s_state << '\n';

	fout << species_vec.size() << '\n';
	for (auto spp : species_vec) spp->save(fout);

	odeStepper.save(fout);
}


void Solver::restore(std::istream &fin, vector<Species_Base*> spp_proto){
	string s; fin >> s; // version number (discard)
	assert(s == "Solver::v2");
	int m;
	fin >> m
	    >> n_statevars_system
		>> current_time
		>> use_log_densities
		>> N0;
	fin >> pi0;
	method = PSPM_SolverType(m);

	fin >> s_state;

	int nspecies;
	fin >> nspecies;
	assert(nspecies == spp_proto.size());
	species_vec = spp_proto;

	for (auto spp : species_vec){
		spp->restore(fin);
	}
	resizeStateFromSpecies(); // this includes system vars
	copyCohortsToState(); // this will overwrite species state (except system vars), but I guess this is better for sake of consistency
	
	odeStepper.restore(fin);
}


// void Solver::printCohortVector(){

// 	std::ofstream cohortprint;
// 	cohortprint.open(std::string("cohort_vector.txt").c_str());

// 	int i = 0;
// 	for (auto s : species_vec){
// 		s->printCohortVector(i, current_time, cohortprint);		
// 		++i;
// 	}

// 	cohortprint.close();
// }

// void Solver::printCohortVector(std::ostream &out){
// 	int i = 0;
// 	for (auto s : species_vec){
// 		s->printCohortVector(i, current_time, out);		
// 		++i;
// 	}
// }



// void Solver::printODEmethod(){
// 	odeStepper.printODEsolvermethod();
// }
