#include <cassert>
#include <algorithm>
#include <numeric>
#include <iomanip>

// *************** Species_Base   ***************  

template<class T>
void Species_Base::addCohort(T bc){
	(dynamic_cast<Species<T>>(*this)).addCohort(bc);
}


// *************** Species<Model> ******************

// template<class Model>
// Species<Model>* Species<Model>::create(){
// 	return new Species<Model>();
// }


template<class Model>
void Species<Model>::resize(int _J){
	J = _J;
	std::cout << "In resize before resize " <<std::endl;
	cohorts.resize(J, boundaryCohort);  // when resizing always insert copies of boundary cohort
	std::cout << "In resize after resize " <<std::endl;
}

template<class Model>
double Species<Model>::get_maxSize(){ // TODO ALERT: make sure this sees the latest state
	if (!X.empty()) return *x.rbegin();	// for FMU, get this from X - TODO: will handle FMU later but right no leaving as is. still cohorts.empty() should probably be first argument
	else if (cohorts.empty()) return xb;
	else {								// else get from state vector
		return cohorts[0].x;			// cohorts are sorted descending
	}
}


// TODO: make this comment more sensible
// Returns an "artificial" vector of largest xk in all the cohorts - no cohort has to actually occupy this point though 

template<class Model>
std::vector<double> Species<Model>::get_maxSizeN(){ // TODO ALERT: make sure this sees the latest state
	if (cohorts.empty()) {
		return xnb;
	}
	else {								// else get from state vector
		std::vector<double> largest = cohorts[0].xn;
		for (size_t i = 1; i <= J; ++i){ //should I check the boundary cohort too? probably not needed can be i < J
			for(size_t k = 0; k < largest.size(); ++k){
					if(cohorts[0].xn[k] > largest[k]){
						largest[k] = cohorts[0].xn[k];
					}
			}
		}
		return largest;			// cohorts are sorted descending
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
	boundaryCohort = Cohort<Model>(M); 
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
	print_extra();
	//if (!X.empty()){
	//    iset.push_back("_X", X.begin(),1);
	//    iset.push_back("_h", h.begin(),1);
	std::cout << "x (" << x.size() << "): "; 
	for (auto xx : x) std::cout << xx << " ";
	std::cout << "\n";
	//}

	std::cout << "Cohorts: (" << cohorts.size() << ")\n";
	std::cout << std::setw(6) << "t0" 
	          << std::setw(12) << "x" 
			  << std::setw(12) << "u" << " ";
	if (!cohorts.empty()){
		for (auto s : cohorts[0].varnames) std::cout << std::setw(12) << s << " ";
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
auto Species<Model>::getXn(int i){
	if (i == -1) return boundaryCohort.xn;
	else return cohorts[i].xn;
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
	// boundaryCohort.set_size(_xb);
}

template <class Model>
void Species<Model>::set_xnb(std::vector<double> _xnb){
	xnb = _xnb;
	boundaryCohort.xn = _xnb;
	boundaryCohort.set_size(_xnb);
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
void Species<Model>::setXn(int i, std::vector<double> _xn){
	cohorts[i].xn = _xn;
	cohorts[i].set_size(_xn);
}

template <class Model>
void Species<Model>::setU(int i, double _u){
	cohorts[i].u = _u;
}

template <class Model>
std::vector <double> Species<Model>::getStateAt(int i){
	std::vector<double> _xn = cohorts[i].xn;
	return(_xn);
}

//To remove
template <class Model>
double Species<Model>::dXn(int i){
	double dxn = 1;
	for (size_t k = 0; k < xn[i].size(); ++k){
		double dx_k = (xn[i][k] - next_xk_desc(xn[i][k], k));
		dxn = dxn * dx_k;
	}
	return dxn;
}

template <class Model> 
double Species<Model>::dXn(std::vector<double> xn1, std::vector<double> xn2){
	double _dxn = 1;
	for (size_t k = 0; k < xn1.size(); ++k){
		double dx_k = (xn2[k] - xn1[k]);
		_dxn = _dxn * dx_k;
	}
	return _dxn;
}

template <class Model> 
std::vector<double> Species<Model>::cohort_dist(std::vector<double> xn1, std::vector<double> xn2){
	std::vector<double> _dxn;
	for (size_t k = 0; k < xn1.size(); ++k){
		double dx_k = (xn2[k] - xn1[k]);
		_dxn.push_back(dx_k);
	}
	return _dxn;
}

template <class Model>
double Species<Model>::next_xk_desc(double xnk, int k){
	double next_smallest = xnb[k];
	for(size_t i = 0; i < xn.size(); i++){
		if(xn[i][k] < xnk && xn[i][k] > next_smallest){
			next_smallest = xn[i][k];
		}
	}
	return next_smallest;
}

template <class Model>
double Species<Model>::next_xk_asc(double xnk, int k){
	double next_biggest = get_maxSizeN()[k];
	for(size_t i = 0; i < xn.size(); i++){
		if(xn[i][k] > xnk && xn[i][k] < next_biggest){
			next_biggest = xn[i][k];
		}
	}
	return next_biggest;
}

// FIXME: need to rethink this... 
template <class Model>
std::vector<double> Species<Model>::next_xn_desc(std::vector<double> _xni){
	std::vector<double> next_smallest = xnb; //boundary is always smaller
	for (size_t i = 1; i <= J; ++i){ //should I check the boundary cohort too? probably not needed can be i < J
		bool smaller = true;
		for(size_t k = 0; k < next_smallest.size(); ++k){
				if(xn[i][k] < _xni[k] && xn[i][k] > next_smallest[k]){
					smaller = false; 
					break;
				}
		}
		if(smaller){
			next_smallest = cohorts[i].xn;
		}
	}
	return next_smallest;
}

// FIXME: need to rethink this... 
template <class Model>
std::vector<double> Species<Model>::next_xn_asc(std::vector<double> _xni){
	std::vector<double> next_biggest = get_maxSizeN();
	for (size_t i = 1; i <= J; ++i){ //should I check the boundary cohort too? probably not needed can be i < J
		bool bigger = true;
		for(size_t k = 0; k < next_biggest.size(); ++k){
			if(xn[i][k] > _xni[k] && xn[i][k] < next_biggest[k]){
				bigger = false; 
				break;
			}
		}
		if(bigger){
			next_biggest = cohorts[i].xn;
		}
	}
	return next_biggest;
}

template <class Model>
void Species<Model>::initExtraState(double t, void * env){
	// init boundary cohort 
	boundaryCohort.init_state(t, env); 
	// init internal cohorts
	for (auto& c : cohorts){
		c.init_state(t, env);	// init state
	}
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
	// std::cout << " ncoh = " << cohorts.size() << "/" << i << std::endl; 
	return cohorts[i].init_density(_x, _env, birth_flux_in);
}

template <class Model>
double Species<Model>::init_density(int i, std::vector<double> _xn, void * env){
	return cohorts[i].init_density(_xn, env, birth_flux_in);
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
std::vector<double> Species<Model>::growthRate(int i, std::vector<double> x, double t, void * env){
	Cohort<Model> &c = (i<0)? boundaryCohort : cohorts[i];
	return c.growthRate(c.xn,t,env);
}


template <class Model>
double Species<Model>::growthRateOffset(int i, double x, double t, void * env){
	Cohort<Model> coff = (i<0)? boundaryCohort : cohorts[i];
	coff.set_size(x);
	//coff.preCompute(coff.x,t,env);

	return coff.growthRate(coff.x,t,env);
}

template <class Model>
double Species<Model>::growthRateOffset(int i, std::vector<double> x, double t, void * env){
	Cohort<Model> coff = (i<0)? boundaryCohort : cohorts[i];
	coff.set_size(x);
	//coff.preCompute(coff.x,t,env);

	return coff.growthRate(coff.xn,t,env);
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
std::vector<double> Species<Model>::growthRateGradient(int i, std::vector<double> x, double t, void * env, std::vector<double> grad_dx){

	// std::cout << "Species::growthRateGradient" << std::endl;
	// std::cout << "x:\t" << x << std::endl;
	// std::cout << "i:\t" << i << std::endl;
	// std::cout << "grad_dx:\t" << grad_dx << std::endl;
	Cohort<Model> &c = (i<0)? boundaryCohort : cohorts[i];

	// std::cout << "c.xn:\t" << c.xn << std::endl;


	std::vector<double> g = c.growthRate(c.xn,t,env);
	// std::cout << "g = c.growthRate(c.xn,t,env):\t" << g << std::endl; 



	std::vector<double> gplus;

	for(int k= 0; k < x.size(); ++k){
		Cohort<Model> cplus = c;
		std::vector<double> _x = x;
		// std::cout << "k:\t" << k << std::endl;
		// std::cout << "_x:\t" << _x << std::endl;
		// std::cout << "grad_dx[k]:\t" << grad_dx[k] << std::endl;

		_x[k] = _x[k] + grad_dx[k];
		// std::cout << "_x:\t" << _x << std::endl;
		cplus.set_size(_x);

		// std::cout << "cplus.xn:\t" << cplus.xn << std::endl;
		std::vector<double> gplus_k = cplus.growthRate(cplus.xn, t, env);

		// std::cout << "gplus_k:\t" << gplus_k[k] << std::endl;
		gplus.push_back((gplus_k[k]-g[k])/grad_dx[k]);
	}
	
	std::vector<double> out;
	// out.insert(g.begin(), g.end());
	out.insert(out.end(), gplus.begin(), gplus.end());
	return out;
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
double Species<Model>::mortalityRate(int i, std::vector<double> x, double t, void * env){
	assert(i>=0);
	Cohort<Model> &c = cohorts[i];
	return c.mortalityRate(c.xn,t,env);
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
std::vector<double> Species<Model>::mortalityRateGradient(int i, std::vector<double> x, double t, void * env, std::vector<double> grad_dx){
	Cohort<Model> &c = (i<0)? boundaryCohort : cohorts[i];

	double g = c.mortalityRate(c.xn,t,env);
	std::vector<double> gplus;

	for(int k= 0; k < x.size(); ++k){
		Cohort<Model> cplus = c;
		std::vector<double> _x = x;
		_x[k] += grad_dx[k];
		cplus.set_size(_x);
		double gplus_k = cplus.mortalityRate(cplus.xn, t, env);
		gplus.push_back((gplus_k-g)/grad_dx[k]);
	}
	
	std::vector<double> out;
	out.push_back(g);
	out.insert(out.end(), gplus.begin(), gplus.end());
	return out;
}
	

template <class Model>
double Species<Model>::birthRate(int i, double x, double t, void * env){
	Cohort<Model> &c = (i<0)? boundaryCohort : cohorts[i];
	return c.birthRate(c.x,t,env);
}

template <class Model>
double Species<Model>::birthRate(int i, std::vector<double> _xn, double t, void * env){
	Cohort<Model> &c = (i<0)? boundaryCohort : cohorts[i];
	return c.birthRate(c.xn,t,env);
}	

template <class Model>
void Species<Model>::getExtraRates(std::vector<double>::iterator &it){
	for (auto& c : cohorts) c.get_rates(it);
}


template <class Model>
void Species<Model>::addCohort(int n){
	// std::cout << "In add cohort, adding n " << n << std::endl;

	if (n > cohorts.max_size()){
		std::cout << "requested n = " << n << " is greater than max_size = " << cohorts.max_size() << std::endl;
	}

	// std::cout << "In add cohort before going through n" << 	std::endl;

	// std::cout << xn << std::endl;
	// std::cout << J << std::endl;

	cohorts.reserve(cohorts.size()+n);
	for (int i=0; i<n; ++i){
		cohorts.push_back(boundaryCohort);
		++J;
	}
	
	// std::cout << "In add cohort after going through n" << 	std::endl;

	// std::cout << cohorts.size() << std::endl;
	// std::cout << J << std::endl;
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

	//std::cout << "Removing cohort no. " << i_min << std::endl;
	cohorts.erase(cohorts.begin()+i_min);
	--J;
}

template <class Model>
void Species<Model>::removeDensestCohortN(){
	if (cohorts.size() < 3) return; // do nothing if there are 2 or fewer cohorts
	int i_min = 0;
	double dx_min = dXn(next_xn_asc(cohorts[i_min].xn), next_xn_desc(cohorts[i_min].xn));
	std::vector<double> maxCohort = get_maxSizeN();
	if(cohorts[i_min].xn == maxCohort){ //this should be executed at most twice
		i_min = ++i_min;
		double dx_min = dXn(next_xn_asc(cohorts[i_min].xn), next_xn_desc(cohorts[i_min].xn));
	}

	for (int i=(i_min + 1); i<J-1; ++i){ // skip first and last cohorts
		if(maxCohort == cohorts[i].xn){
			continue;
		}
		double dx = dXn(next_xn_asc(cohorts[i].xn), next_xn_desc(cohorts[i].xn));
		if (dx < dx_min){
			dx_min = dx;
			i_min = i;
		}
	}

	//std::cout << "Removing cohort no. " << i_min << std::endl;
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
void Species<Model>::removeDenseCohorts(std::vector<double> dxcut){
	// mark cohorts to remove; skip 1st and last cohort
	if (cohorts.size() < 3) return; // do nothing if there are 2 or fewer cohorts
	std::vector<double> maxCohort = get_maxSizeN();
	
	for (int i=0; i<J-1; i+=2){
		if(maxCohort == cohorts[i].xn){
			continue;
		}
		std::vector<double> dx_lo = cohort_dist(next_xn_asc(cohorts[i].xn), cohorts[i].xn);
		std::vector<double> dx_hi = cohort_dist(cohorts[i].xn, next_xn_desc(cohorts[i].xn));

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
void Species<Model>::mergeCohortsAddU(double dxcut){
	// mark cohorts to remove; skip 1st and last cohort
	for (int i=1; i<J-1; ++i){
		double dx = cohorts[i-1].x-cohorts[i].x;

		if (dx < dxcut){
			cohorts[i-1].remove = true;
			// FIXME: Need to also average extra state?
			cohorts[i].x = (cohorts[i].x*cohorts[i].u + cohorts[i-1].x*cohorts[i-1].u)/(cohorts[i].u + cohorts[i-1].u);
			cohorts[i].u = cohorts[i].u + cohorts[i-1].u;
		}
	}

	// remove marked cohorts
	auto it_end = std::remove_if(cohorts.begin(), cohorts.end(), [](Cohort<Model> &c){return c.remove;});
	cohorts.erase(it_end, cohorts.end());

	// reset size
	J = cohorts.size();

}

template <class Model>
void Species<Model>::mergeCohortsAddU(std::vector<double> dxcut){
	// mark cohorts to remove; skip 1st and last cohort
	std::vector<double> maxCohort = get_maxSizeN();
	for (int i=0; i<J-1; i+=2){
		if(maxCohort == cohorts[i].xn){
			continue;
		}
		std::vector<double> dx = cohort_dist(cohorts[i].xn, next_xn_asc(cohorts[i].xn));
		if (dx < dxcut){
			cohorts[i-1].remove = true;
			// FIXME: Need to also average extra state?
			for(int k=0; k<xnb.size(); ++k){
				cohorts[i].xn[k] = (cohorts[i].xn[k]*cohorts[i].u + cohorts[i-1].xn[k]*cohorts[i-1].u)/(cohorts[i].u + cohorts[i-1].u);
			}

			cohorts[i].u = cohorts[i].u + cohorts[i-1].u;
		}
	}

	// remove marked cohorts
	auto it_end = std::remove_if(cohorts.begin(), cohorts.end(), [](Cohort<Model> &c){return c.remove;});
	cohorts.erase(it_end, cohorts.end());

	// reset size
	J = cohorts.size();

}



template <class Model>
void Species<Model>::sortCohortsDescending(int skip){
	std::sort(cohorts.begin(), cohorts.end()-skip, [](const Cohort<Model> &a, const Cohort<Model> &b){return a.x > b.x;});
}


template <class Model>
void Species<Model>::save(std::ofstream &fout){
//	Species_Base::save(fout);
	fout << "Species<T>::v1\n";
	fout << std::make_tuple(
		J
	  , n_extra_statevars
	  , noff_abm
	  , birth_flux_in
	  , bfin_is_u0in
	  , xb);
	fout << '\n';
	fout << X << x << h;

	boundaryCohort.save(fout, n_extra_statevars);
	for (auto& C : cohorts) C.save(fout, n_extra_statevars);
}

template <class Model>
void Species<Model>::restore(std::ifstream &fin){
//	Species_Base::restore(fin);
	std::cout << "Restoring Species<T>" << std::endl;
	std::string s; fin >> s; // version number (discard)
	assert(s == "Species<T>::v1");
	fin >> J	
	    >> n_extra_statevars
	    >> noff_abm
	    >> birth_flux_in
	    >> bfin_is_u0in
	    >> xb;
	
	fin >> X >> x >> h;

	boundaryCohort.restore(fin, n_extra_statevars);
	std::cout << "In restore before resize " <<std::endl;
	cohorts.resize(J, boundaryCohort); // cohorts must always be copy-constructed from the boundary cohort

	std::cout << "In restore after resize " <<std::endl;
	for (auto& C : cohorts) C.restore(fin, n_extra_statevars);
}


//TODO: maybe fix this later
template <class Model>
void Species<Model>::printCohortVector(int speciesInd, double time, std::ostream &out){

	// std::cout << "J is " << J << std::endl;
	for(int i=0; i < cohorts.size(); ++i){
		cohorts[i].print(speciesInd, time, out);	
		out << "\n";	
	}

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


