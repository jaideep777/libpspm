#include <cassert>
#include <algorithm>
#include <numeric>

// *************** Species_Base   ***************  

template<class T>
void Species_Base::addCohort(T bc){
	(dynamic_cast<Species<T>>(*this)).addCohort(bc);
}


// *************** Species<Model> ******************



template<class Model>
void Species<Model>::resize(int _J){
	J = _J;
	cohorts.resize(J, boundaryCohort);  // when resizing, always insert copies of boundary cohort
}


template<class Model>
double Species<Model>::get_maxSize(){ // TODO ALERT: make sure this sees the latest state
	if (!X.empty()) return *x.rbegin();	// for FMU, get this from X
	else if (cohorts.empty()) return 0;
	else {								// else get from state vector
		return cohorts[0].x;			// cohorts are sorted descending
	}
}


// FIXME: this constructor should  be corrected or removed.
template<class Model>
Species<Model>::Species(std::vector<double> breaks){
	J = breaks.size(); // changes with 
	cohorts.resize(J, boundaryCohort);
	for (int i=0; i<J; ++i) cohorts[i].x = breaks[i];	
}


template<class Model>
Species<Model>::Species(Model M){
	boundaryCohort = savedCohort = Cohort<Model>(M); 
}


template <class Model>
void Species<Model>::print(){
	std::cout << "~~~~~ Species ~~~~~\n";
	//auto iset = get_iterators(sv);
	//std::cout << "start index = " << start_index <<"\n";
	//std::cout << "Model = " << mod << "\n";
	std::cout << "xb = " << boundaryCohort.x << " / " << xb << "\n";
	std::cout << "xsize = " << J << "\n";
	std::cout << "Extra state variables: " << n_extra_statevars << "\n";
	std::cout << "Input birth flux = " << birth_flux_in << "\n";
	//if (!X.empty()){
	//    iset.push_back("_X", X.begin(),1);
	//    iset.push_back("_h", h.begin(),1);
	std::cout << "x (" << x.size() << "): "; 
	for (auto xx : x) std::cout << xx << " ";
	std::cout << "\n";
	//}

	std::cout << "Cohorts: (" << cohorts.size() << ")\n";
	std::cout << "t0" << "\tx" << "\t" << "u" << "\t";
	if (!cohorts.empty()){
		for (auto s : cohorts[0].varnames) std::cout << s << "\t";
	}
	std::cout << "\n";
	for (auto& c : cohorts) {
		//c.print_xu(); //std::cout << c.x << "\t" << c.u << "\t";
		c.print();
		std::cout << "\n";
	}
	std::cout << "- - - - - - - - - - - - - - - - - - - - - - - \n";
	//boundaryCohort.print_xu(); //std::cout << c.x << "\t" << c.u << "\t";
	boundaryCohort.print();
	std::cout << "\n";
	std::cout << "Max size = " << get_maxSize() << "\n";
	//std::cout << "State (" << size() << "):\n";
	//iset.print();
	
	//std::cout << "Rates (" << size() << "):\n";
	//auto irates = get_iterators(rv);
	//irates.print();
	std::cout << "-------\n\n"; std::cout.flush();

}

template <class Model>
Cohort<Model>& Species<Model>::getCohort(int i){
	if (i == -1) return boundaryCohort;
	else return cohorts[i];
}

template <class Model>
double Species<Model>::getX(int i){
	if (i == -1) return boundaryCohort.x;
	else return cohorts[i].x;
}


template <class Model>
double Species<Model>::getU(int i){
	return cohorts[i].u;
}


template <class Model>
double Species<Model>::get_boundary_u(){
	return boundaryCohort.u;
}


template <class Model>
void Species<Model>::set_xb(double _xb){
	xb = _xb;
	boundaryCohort.x = _xb;
	boundaryCohort.set_size(_xb);
}


template <class Model>
void Species<Model>::set_ub(double _ub){
	boundaryCohort.u = _ub;
}


template <class Model>
void Species<Model>::set_birthTime(int i, double t0){
	cohorts[i].birth_time = t0;
}


template <class Model>
void Species<Model>::setX(int i, double _x){
	cohorts[i].x = _x;
	cohorts[i].set_size(_x);
}


template <class Model>
void Species<Model>::setU(int i, double _u){
	cohorts[i].u = _u;
}


template <class Model>
void Species<Model>::initAndCopyExtraState(double t, void * env, std::vector<double>::iterator &it){
	// init boundary cohort (no copy required)
	boundaryCohort.init_state(t, env); 
	// init internal cohorts and copy to state vector
	for (auto& c : cohorts){
		c.init_state(t, env);	// init state
		auto it_prev = it;		
		it = c.get_state(it);	// copy the initialized state into state vector
		assert(distance(it_prev, it) == n_extra_statevars);
	}
}


template <class Model>
void Species<Model>::initBoundaryCohort(double t, void * env){
	boundaryCohort.birth_time = t;
	boundaryCohort.init_state(t, env);
}


template <class Model>
double Species<Model>::init_density(int i, double _x, void * _env){
	return cohorts[i].init_density(_x, _env, birth_flux_in);
}

// TODO: check increment here itself
template <class Model>
void Species<Model>::copyExtraStateToCohorts(std::vector<double>::iterator &it){
	for (auto& c : cohorts) c.set_state(it);
}


template <class Model>
void Species<Model>::copyCohortsExtraToState(std::vector<double>::iterator &it){
	for (auto& c : cohorts) c.get_state(it);
}


template <class Model>
void Species<Model>::triggerPreCompute(){
	for (auto& c : cohorts) c.need_precompute = true;
	boundaryCohort.need_precompute = true;
}


template <class Model>
double Species<Model>::establishmentProbability(double t, void * env){
	return boundaryCohort.establishmentProbability(t, env);
}


// FIXME: Should be renamed as "calc_boundary_u"
template <class Model>
double Species<Model>::calc_boundary_u(double gb, double pe){
	//std::cout << "calc_boundary_u\n";
	if (bfin_is_u0in){
		boundaryCohort.u = birth_flux_in;
	}
	else {
		boundaryCohort.u = (gb>0)? birth_flux_in * pe/gb  :  0; 
	}
	return boundaryCohort.u;
}


template <class Model>
double Species<Model>::growthRate(int i, double x, double t, void * env){
	Cohort<Model> &c = (i<0)? boundaryCohort : cohorts[i];
	return c.growthRate(c.x,t,env);
}


template <class Model>
double Species<Model>::growthRateOffset(int i, double x, double t, void * env){
	Cohort<Model> coff = (i<0)? boundaryCohort : cohorts[i];
	coff.set_size(x);
	//coff.preCompute(coff.x,t,env);

	return coff.growthRate(coff.x,t,env);
}


template <class Model>
std::vector<double> Species<Model>::growthRateGradient(int i, double x, double t, void * env, double grad_dx){
	Cohort<Model> &c = (i<0)? boundaryCohort : cohorts[i];

	Cohort<Model> cplus = c;
	cplus.set_size(c.x + grad_dx);
	//cplus.preCompute(cplus.x,t,env);
	
	double g = c.growthRate(c.x,t,env);
	double gplus = cplus.growthRate(cplus.x, t, env);

	return {g, (gplus-g)/grad_dx};
}


template <class Model>
std::vector<double> Species<Model>::growthRateGradientCentered(int i, double xplus, double xminus, double t, void * env){
	assert(i>=0);

	Cohort<Model> cplus = cohorts[i];
	cplus.set_size(xplus);
	//cplus.preCompute(cplus.x,t,env);

	Cohort<Model> cminus = cohorts[i];
	cminus.set_size(xminus);
	//cminus.preCompute(cminus.x,t,env);
	
	double gplus  = cplus.growthRate(cplus.x, t, env);
	double gminus = cminus.growthRate(cminus.x, t, env);
	
	return {gplus, gminus};
}


template <class Model>
double Species<Model>::mortalityRate(int i, double x, double t, void * env){
	assert(i>=0);
	Cohort<Model> &c = cohorts[i];
	return c.mortalityRate(c.x,t,env);
}
	

template <class Model>
std::vector<double> Species<Model>::mortalityRateGradient(int i, double x, double t, void * env, double grad_dx){
	Cohort<Model> &c = (i<0)? boundaryCohort : cohorts[i];

	Cohort<Model> cplus = c;
	cplus.set_size(x + grad_dx);
	//cplus.preCompute(cplus.x,t,env);
	
	double g = c.mortalityRate(c.x,t,env);
	double gplus = cplus.mortalityRate(cplus.x, t, env);
	
	return {g, (gplus-g)/grad_dx};
}
	

template <class Model>
double Species<Model>::birthRate(int i, double x, double t, void * env){
	Cohort<Model> &c = (i<0)? boundaryCohort : cohorts[i];
	return c.birthRate(c.x,t,env);
}
	

template <class Model>
void Species<Model>::getExtraRates(std::vector<double>::iterator &it){
	for (auto& c : cohorts) c.get_rates(it);
}


template <class Model>
void Species<Model>::addCohort(int n){
	cohorts.reserve(cohorts.size()+n);
	for (int i=0; i<n; ++i){
		cohorts.push_back(boundaryCohort);
		++J;
	}
}



// This function allows a user to add a custom cohort to the species any time.
// However, Solver->resizeStateFromSpecies() must be called after calling this function.
template <class Model>
void Species<Model>::addCohort(Cohort<Model> bc){
	cohorts.push_back(bc);
	++J;
	// FIXME: add option to sort cohorts here.
}


template <class Model>
void Species<Model>::markCohortForRemoval(int i){
	cohorts[i].remove = true;
}

template <class Model>
void Species<Model>::removeMarkedCohorts(){
	// remove marked cohorts
	auto it_end = std::remove_if(cohorts.begin(), cohorts.end(), [](Cohort<Model> &c){return c.remove;});
	cohorts.erase(it_end, cohorts.end());
	
	// reset size
	J = cohorts.size();
}


template <class Model>
void Species<Model>::removeDensestCohort(){
	if (cohorts.size() < 3) return; // do nothing if there are 2 or fewer cohorts
	int i_min = 1;
	double dx_min = cohorts[0].x - cohorts[2].x;  
	for (int i=1; i<J-1; ++i){ // skip first and last cohorts
		double dx = cohorts[i-1].x - cohorts[i+1].x;
		if (dx < dx_min){
			dx_min = dx;
			i_min = i;
		}
	}

	cohorts.erase(cohorts.begin()+i_min);
	--J;
}


template <class Model>
void Species<Model>::removeDenseCohorts(double dxcut){
	// mark cohorts to remove; skip 1st and last cohort
	for (int i=1; i<J-1; i+=2){
		double dx_lo = cohorts[i-1].x-cohorts[i].x;
		double dx_hi = cohorts[i].x-cohorts[i+1].x;

		if (dx_lo < dxcut || dx_hi < dxcut) cohorts[i].remove = true;
	}

	// remove marked cohorts
	auto it_end = std::remove_if(cohorts.begin(), cohorts.end(), [](Cohort<Model> &c){return c.remove;});
	cohorts.erase(it_end, cohorts.end());

	// reset size
	J = cohorts.size();
}


template <class Model>
void Species<Model>::removeDeadCohorts(double ucut){
	// mark cohorts to remove; skip pi0-cohort (index J-1)
	for (int i=0; i<J-1; ++i){
		if (cohorts[i].u < ucut) cohorts[i].remove = true;
	}

	// remove marked cohorts
	auto it_end = std::remove_if(cohorts.begin(), cohorts.end(), [](Cohort<Model> &c){return c.remove;});
	cohorts.erase(it_end, cohorts.end());

	// reset size
	J = cohorts.size();
}


template <class Model>
void Species<Model>::sortCohortsDescending(){
	std::sort(cohorts.begin(), cohorts.end(), [](const Cohort<Model> &a, const Cohort<Model> &b){return a.x > b.x;});
}

//template <class Model>
//void Species<Model>::backupCohort(int j){
//	savedCohort = cohorts[j];
//}

//template <class Model>
//void Species<Model>::restoreCohort(int j){
//	cohorts[j] = savedCohort;
//}

//template <class Model>
//void Species<Model>::copyBoundaryCohortTo(int j){
//	cohorts[j] = boundaryCohort;
//}

//template <class Model>
//void Species<Model>::backupBoundaryCohort(){
	//boundaryCohort_backup = boundaryCohort;
//}

//template <class Model>
//void Species<Model>::restoreBoundaryCohort(){
	//boundaryCohort = boundaryCohort_backup;
//}


