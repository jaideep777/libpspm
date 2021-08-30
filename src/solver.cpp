#include "solver.h"

#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include <algorithm>
#include <functional>

#include "iterator_set.h"

using namespace std;


std::vector <double> seq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

std::vector <double> logseq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = exp(log(from) + i*(log(to)-log(from))/(len-1));
	return x;
}
// ~~~~~~~~~~~ SOLVER ~~~~~~~~~~~~~~~~~~~~~

//template<class Model, class Environment>
//const int Solver<Model,Environment>::xsize(){
//    if (method == SOLVER_FMU) return J;	
//    if (method == SOLVER_MMU) return J;  
//    if (method == SOLVER_CM ) return J;
//    if (method == SOLVER_EBT) return J;
//}

Solver::Solver(PSPM_SolverType _method) : odeStepper(0, 1e-6, 1e-6) {
	method = _method;
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

	if (method == SOLVER_FMU){
		std::sort(xbreaks.begin(), xbreaks.end(), std::less<double>());  // sort cohorts ascending for FMU
		s->xb = *xbreaks.begin();
	}
	else {
		std::sort(xbreaks.begin(), xbreaks.end(), std::greater<double>());  // sort cohorts descending for CM, EBT
		s->xb = *xbreaks.rbegin();
	}

	int J;
	if (method == SOLVER_FMU) J = xbreaks.size()-1;	
	if (method == SOLVER_MMU) J = xbreaks.size()-1;  
	if (method == SOLVER_CM ) J = xbreaks.size();
	if (method == SOLVER_EBT) J = xbreaks.size();

	s->x = xbreaks;
	s->resize(J);

	if (method == SOLVER_FMU) n_statevars_internal = 1;
	else n_statevars_internal = 2;

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


void Solver::resetState(){  // FIXME: This is currently redundant, and needs to be improved with reset of both state and cohorts for a true reset of state
	current_time = 0;
	odeStepper = RKCK45<vector<double>> (0, control.ode_eps, control.ode_initial_step_size);  // this is a cheap operation, but this will empty the internal containers, which will then be (automatically) resized at next 1st ODE step. Maybe add a reset function to the ODE stepper? 

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
	string types[] = {"FMU", "MMU", "CM", "EBT"};
	std::cout << "+ Type: " << types[method] << std::endl;

	std::cout << "+ State size = " << state.size() << "\n";
	std::cout << "+ Rates size = " << rates.size() << "\n";
	std::cout << "+ Species:\n";
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
		if (method == SOLVER_FMU){
			for (size_t i=0; i<s->J; ++i){
				double X = (s->x[i]+s->x[i+1])/2;			// for FMU, X is the midpoint of the cell edges
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
		
		if (method == SOLVER_FMU){
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
		
		if (method == SOLVER_FMU){
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


// current_time is updated by the ODE solver at every (internal) step
void Solver::step_to(double tstop){
	// do nothing if tstop is <= current_time
	if (tstop <= current_time) return;
	
	if (method == SOLVER_FMU){	
		auto derivs = [this](double t, vector<double> &S, vector<double> &dSdt){
			copyStateToCohorts(S.begin());
			env->computeEnv(t, this);
			// precompute all species (prepare for rate calcs)
			for (int k = 0; k<species_vec.size(); ++k) preComputeSpecies(k,t);	

			this->calcRates_FMU(t, S, dSdt);
		};
		
		odeStepper.Step_to(tstop, current_time, state, derivs); // state = [U]
		copyStateToCohorts(state.begin());
	}
	
	
	//if (method == SOLVER_MMU){
	//}
	
	
	if (method == SOLVER_EBT){
		auto derivs = [this](double t, vector<double> &S, vector<double> &dSdt){
			// copy state vector to cohorts
			copyStateToCohorts(S.begin());
			
			// compute environment
			env->computeEnv(t, this);

			// precompute all species (prepare for rate calcs)
			for (int k = 0; k<species_vec.size(); ++k) preComputeSpecies(k,t);	

			// get rates
			calcRates_EBT(t, S, dSdt);
		};
		
		// integrate 
		odeStepper.Step_to(tstop, current_time, state, derivs); // state = [pi0, Xint, N0, Nint]
		
		// after the last ODE step, the state vector is updated but cohorts still hold an intenal ODE state (y+k5*h etc).
		// normally, this will be no problem since state will be copied to cohorts in the next rates call. 
		// But since add/remove cohort below will rewrite the state from cohorts, the updated state vector will be lost
		// rewrite the cohorts now to avoid this.
		copyStateToCohorts(state.begin());
		
		// update cohorts
		removeDeadCohorts_EBT();
		addCohort_EBT();  // Add new cohort if N0 > 0. Add after removing dead ones otherwise this will also be removed. 
	}
	
	
	if (method == SOLVER_CM){
		auto derivs = [this](double t, vector<double> &S, vector<double> &dSdt){
			// copy state vector to cohorts
			copyStateToCohorts(S.begin());

			// update u0 (u of boundary cohort)
			for (auto s : species_vec) s->get_u0(t, env);
			
			// compute environment
			env->computeEnv(t, this);

			// precompute all species (prepare for rate calcs)
			for (int k = 0; k<species_vec.size(); ++k) preComputeSpecies(k,t);	

			// get rates
			calcRates_CM(t, S, dSdt);
		};
		
		// integrate 
		odeStepper.Step_to(tstop, current_time, state, derivs); // state = [pi0, Xint, N0, Nint]
		
		// after the last ODE step, the state vector is updated but cohorts still hold an intenal ODE state (y+k5*h etc).
		// normally, this will be no problem since state will be copied to cohorts in the next rates call. 
		// But since add/remove cohort below will rewrite the state from cohorts, the updated state vector will be lost
		// rewrite the cohorts now to avoid this.
		copyStateToCohorts(state.begin());

		// update cohorts
		if (control.update_cohorts){
			addCohort_CM();		// add before so that it becomes boundary cohort and first internal cohort can be (potentially) removed
			removeCohort_CM();
		}
		//env->computeEnv(current_time, this); // is required here IF rescaleEnv is used in derivs
	}
}


void Solver::preComputeSpecies(int k, double t){
	auto spp = species_vec[k];

	double pi0, N0;
	// backup and real-ize pi0-cohort
	if (method == SOLVER_EBT){ // for EBT, we need to pi0-cohort too.
		// get pi0, N0 from last cohort
		pi0  =  spp->getX(spp->J-1);
		N0   =  spp->getU(spp->J-1);
		
		// update pi0-cohort with actual x0 value
		double x0 = spp->xb + pi0/(N0+1e-12);
		spp->setX(spp->J-1, x0);
	}
	
	// precompute cohorts
	spp->preComputeAllCohorts(t,env);

	// restore pi0-cohort
	if (method == SOLVER_EBT){ // for EBT, we need to pi0-cohort too.
		spp->setX(spp->J-1, pi0);
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

	
vector<double> Solver::newborns_out(){
	// update Environment from latest state
	//copyStateToCohorts(state.begin());
	env->computeEnv(current_time, this);
	for (int k = 0; k<species_vec.size(); ++k) preComputeSpecies(k,current_time);	

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
		//double birthFlux = integrate_x(newborns_production, current_time, k);
		//double birthFlux = integrate_wudx_above(newborns_production, current_time, 0, k);
		double birthFlux = calcSpeciesBirthFlux(k, current_time);
		b_out.push_back(birthFlux);
	}
	return b_out;
}

// FOR DEBUG ONLY, using TESTMODEL
vector<double> Solver::u0_out(){
	vector <double> u0out;
	vector <double> newbornsout = newborns_out();
	for (int k=0; k < species_vec.size(); ++k){
		//species_vec[k]->preCompute(-1, current_time, env); // not req because precomputeAllCohorts called in newborns_out() precomputes BC too.	
		u0out.push_back(newbornsout[k]/species_vec[k]->growthRate(-1, species_vec[k]->xb, current_time, env));
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


