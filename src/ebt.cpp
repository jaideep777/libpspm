#include <algorithm>
#include <cassert>

#include <solver.h>

void Solver::calcRates_EBT(double t, vector<double>&S, vector<double> &dSdt){

	vector<double>::iterator its = S.begin()    + n_statevars_system; // Skip system variables
	vector<double>::iterator itr = dSdt.begin() + n_statevars_system;
	
	for (int s = 0; s < species_vec.size(); ++s){
		Species_Base * spp = species_vec[s];

		double   pi0  =  spp->getX(spp->J-1);	 // last cohort is pi0, N0
		double   N0   =  spp->getU(spp->J-1);
		//std::cout << "pi = " << pi0 << ", N0 = " << N0 << "\n";

		double grad_dx = 0.001;
		std::vector<double> g_gx = spp->growthRateGradient(-1, spp->xb, t, env, grad_dx);
		std::vector<double> m_mx = spp->mortalityRateGradient(-1, spp->xb, t, env, grad_dx);
		//std::cout << "g = " << g_gx[0] << ", gx = " << g_gx[1] << "\n";
		//std::cout << "m = " << m_mx[0] << ", mx = " << m_mx[1] << "\n";

		double birthFlux;
		if (spp->birth_flux_in < 0){	
			birthFlux = integrate_x([this, s](int i, double t){
												return species_vec[s]->birthRate(i,species_vec[s]->getX(i),t, env);
											}, t, s);
		}
		else{
			double u0 = spp->get_u0(t, env);
			birthFlux = u0*g_gx[0];
			//double g = spp.mod->growthRate(spp.xb, current_time, env);
			//// --- debug ---
			//if (spp.bfin_is_u0in)
				//u[0] = spp.birth_flux_in;
			//else{
			//// -------------
				//double d = (g>0)? spp.birth_flux_in * spp.mod->establishmentProbability(current_time, env)/g  : 0;
				//u[0] = d; 
			//}
		}
		//cout << "bf = " << birthFlux << "\n";
		//if (spp.birth_flux_in < 0){
			//birthFlux = integrate_x([spp, this](double z, double t){return spp.mod->birthRate(z,t,env);}, t, S, species_id);
		//}
		//else {
			//// --- debug ---
			//if (spp.bfin_is_u0in)
				//birthFlux = (gb>0)? spp.birth_flux_in*spp.mod->establishmentProbability(t, env) * gb : 0;
			//else{
			//// -------------
				//birthFlux = (gb>0)? spp.birth_flux_in*spp.mod->establishmentProbability(t, env)  : 0;
			//}
		//}


		for (int i=0; i<spp->J-1; ++i){	// go down to the second last cohort (exclude boundary cohort)
			*itr++ =  spp->growthRate(i, spp->getX(i), t, env);							// dx/dt
			*itr++ = -spp->mortalityRate(i, spp->getX(i), t, env) * spp->getU(i);		// du/dt
			its += 2;
		}

		// dpi0/dt and dN0/dt
		double mb = m_mx[0], mortGrad = m_mx[1], gb = g_gx[0], growthGrad = g_gx[1];	
		double dN0 = -mb*N0 - mortGrad*pi0 + birthFlux;
		double dpi0  = gb*N0 + growthGrad*pi0 - mb*pi0;
		//std::cout << "dpi0 = " << dpi0 << ", dN0 = " << dN0 << "\n";
		*itr++ = dpi0;
		*itr++ = dN0;
		its += 2;

		if (spp->n_extra_statevars > 0){
			auto itr_prev = itr;
			spp->getExtraRates(itr);
			assert(distance(itr_prev, itr) == spp->n_extra_statevars*spp->J); 
			its += spp->n_extra_statevars*spp->J; 	
		}
		
		// TODO: Should other rates of boundary cohort be explicitly set to zero? YES, set to nan.

		//if (spp.varnames_extra.size() > 0){
		//    auto it_returned = spp.mod->calcRates_extra(*itx, t, env, itse, itre);
		//    assert(distance(itre, it_returned) == spp.varnames_extra.size());
		//}
			
		
	}

}



//template <class Model, class Environment>
void Solver::addCohort_EBT(){
	
	for (auto spp : species_vec){
		// 1. internalize the pi0-cohort
		//   get pi0, N0 from last cohort
		double   pi0  =  spp->getX(spp->J-1);
		double   N0   =  spp->getU(spp->J-1);

		//   replace pi0-cohort with boundary cohort
		spp->copyBoundaryCohortTo(spp->J-1);

		//   update (new) pi0-cohort with actual x0 value
		double x0 = spp->xb + pi0/(N0+1e-12);
		spp->setX(spp->J-1, x0);
		spp->setU(spp->J-1, N0);

		// 2. insert a new cohort
		spp->addCohort(); // this will add a copy of boundary cohort

		// 3. set the new cohort as the pi0-cohort
		spp->setX(spp->J-1, 0); // set pi0 and N0 to 0
		spp->setU(spp->J-1, 0);
	}

	// 4. reset state from cohorts
	resizeStateFromSpecies();
	copyCohortsToState();
	
	//// Create index and value vectors to insert into state
	//vector <int> at;
	//vector <double> vals;
	//int n_insertions = 0;	// keep track of cumulative number of insertions so that species start_index can be updated

	//for (auto& spp : species_vec){
	//    // 1. internalize boundary cohort
	//    auto iset = spp.get_iterators(state);
	//    auto itx = iset.get("X");
	//    auto itu = iset.get("u");
	//    auto itse = iset.get((spp.varnames_extra.size()>0)? spp.varnames_extra[0] : "X");	// dummy init //FIXME: commented line below does not work. Explore why
	
	//    if (*itu == 0) continue;  // skip addition of new cohort if N0 is 0.

	//    // 1a. set size of internalized cohort
	//    *itx = spp.xb + *itx/(*itu+1e-12);
		
	//    // 1b. init extra state variables of internalized cohort
	//    if (spp.varnames_extra.size() > 0){
	//        int id = iset.getIndex(spp.varnames_extra[0]); // get index of 1st extra variable
	//        for (iset.begin(); !iset.end(); ++iset){
	//            vector <double> v = spp.mod->initStateExtra(*itx, current_time, env);  // returned vector will be `move`d so this is fast 
	//            auto it_vec = iset.get();
	//            for (size_t i = 0; i<spp.varnames_extra.size(); ++i){
	//                *it_vec[id+i] = v[i];
	//            }
	//        }
	//    }
		
	//    // 2. add new cohort
	//    int spp_start_index_new = spp.start_index + n_insertions; // calculate new start_index for species and store it

	//    at.push_back(spp.start_index + 0);		// TODO: should these be generalized to use iteratorset?
	//    at.push_back(spp.start_index + spp.J);
	//    vals.push_back(0);
	//    vals.push_back(0);
	//    n_insertions += 2;
	//    for (int i=0; i<spp.varnames_extra.size(); ++i){
	//        at.push_back(spp.start_index + 2*spp.J);
	//        vals.push_back(0);
	//        ++n_insertions;	
	//    }
		
	//    spp.start_index = spp_start_index_new;	// update start_index of species
	//    ++spp.J;								// increment J to reflect new system size
	//    setupLayout(spp);						// rebuild parameters for IteratorSet
	//}

	//vector_insert(state, at, vals);		// insert new cohorts into state
	//rates.resize(state.size(), -999);	// resize rates to correct size

}


//template<class Model, class Environment>
void Solver::removeDeadCohorts_EBT(){

	for (auto spp : species_vec){
		spp->removeDeadCohorts(1e-10);	
	}

	resizeStateFromSpecies();
	copyCohortsToState();

	//// TODO: this should be made truely generic - maybe after making of indexset
	////// cohorts are x0, x1, ...xJ-1,  u0, u1, .... uJ-1, a0,b0,c0,..., a1,b1,c1,...,aJ-1,bJ-1,cJ-1
	////auto px = state.begin(); advance(px, 1); // point at x1
	////auto pu = state.begin(); advance(pu, xsize()+1); // point at u1
	////auto last = state.begin(); advance(last, xsize()-1); // point at xJ (1 past the last value to be considered)

	//////cout << *px << " " << *last << " " << *pu << endl;
	//vector<double> flags(state.size(), 0.0);

	//int n_removals = 0;
	//for (auto& spp : species_vec){
	//    //if (spp.J < 500) continue;   // remove a cohort if number of cohorts in the species exceeds a threshold

	//    int spp_start_index_new = spp.start_index - n_removals;  // start_index must be shifted left by number of elements removed so far
	//    auto iset = spp.get_iterators(state);
	//    auto iset_f = spp.get_iterators(flags);
		
	//    auto& itu = iset.get("u");
		
	//    iset.begin(); iset_f.begin();
	//    ++iset; ++iset_f; // point to second element (skip boundary cohort)

	//    //int n_cohorts_removed = 0;
	//    for (; iset.dist < iset.size-1; ++iset, ++iset_f){	// go till second-last element
	//        if (*itu < 1e-6){
	//            --spp.J; //++n_cohorts_removed;
	//            // mark flags to remove
	//            auto iset_vec_f = iset_f.get();  // vector of iterators to remove. FIXME: How to ensure that extra variables dont come up in this vector which are not part of state?
	//            for (int i=0; i<iset_vec_f.size(); ++i){
	//                *iset_vec_f[i] = 1.0;
	//                ++n_removals;
	//            }
	//        }
	//    }
	
	//    spp.start_index = spp_start_index_new;
	//    //spp.J -= n_cohorts_removed; // reduce species size
	//    setupLayout(spp);	
	//}
	
	//// remove all flagged elements from state
	//auto pred = [this, &flags](const double& s) -> bool {
	//    return flags[&s - (const double*)&this->state[0]] == 1.0; 
	//};

	//auto pend = std::remove_if(state.begin(), state.end(), pred);
	//state.erase(pend, state.end());
	//////for (auto z : state) cout << z << " "; cout << endl;

}


//template <class Model>
//void Solver<Model>::addCohort_EBT(){
//    // internalize boundary cohort
//    auto p_pi0  = state.begin();
//    auto p_N0   = state.begin() + xsize(); 
	
//    double x0 = xb + *p_pi0/(*p_N0+1e-12); 
//    *p_pi0 = x0; // pi0 becomes x0 when boundary cohort is internalized
	
//    // insert new boundary cohort
//    vector <int> at = {0, xsize()};
//    vector <double> vals = {0, 0};
//    vector <double> ex = mod->initStateExtra(xb, current_time);
//    for (int i=0; i<varnames_extra.size(); ++i){
//        at.push_back(2*xsize());
//        vals.push_back(ex[i]); 
//    }

//    vector_insert(state, at, vals);
//    rates.resize(state.size(), -999);
//    ++J;	// increment J to reflect new system size
//    setupLayout();  // rebuild parameters for IteratorSet

//}


//// this function was tested here: cpp.sh/67t3i
//template <class Model>
//void Solver<Model>::removeDeadCohorts_EBT(){
//    // since erase() invalidates all iterators beyond erased position, 
//    // we flag all elements to be removed in 1 pass, 
//    // then in a 2nd pass, iterate backwards removing the flagged cohorts 
//    // One could have simply set the values of flagged x and N to some random value (e.g. 1e-20) and use it in remove_if(); but that is not fully failsafe
//    // using iterators so that this code works even if state is a list
//    vector<bool> flags(state.size(), false);
	
//    // iterators to beginning of x and beginning of N in state array
//    auto p_N = state.begin() + xsize();
//    auto p_x_flag = flags.begin();
//    auto p_N_flag = flags.begin() + xsize();

//    // for dead cohorts, mark N and correponding x for removal  
//    ++p_N; ++p_x_flag; ++p_N_flag; // skip pi0, N0 from removal
//    int Jnew = J;
//    for (int i=0; i<xsize()-1; ++i){
//        if (*p_N < 1e-10){
//            *p_x_flag = *p_N_flag = true; 
//            --Jnew;
//        }
//        ++p_N; ++p_x_flag; ++p_N_flag; 
//    }
//    J = Jnew;

//    // remove flagged elements
//    auto pred = [this, &flags](const double& s) -> bool {
//        return flags[&s - (const double*)&this->state[0]]; // TODO: make conatiner-type safe
//    };

//    auto pend = std::remove_if(state.begin(), state.end(), pred);
//    state.erase(pend, state.end());

//}



//// FIXME: Need to be fixed. Cross check with R
//template <class Model>
//vector<double> Solver<Model>::cohortsToDensity_EBT(vector <double> &breaks){

//    auto pX = state.begin();
//    auto pN = state.begin()+xsize();

//    *pX = xb + *pX/(*pN+1e-12); 
				
//    vector<double> dens;
//    vector<double> xmean(1, xb);

//    for (int i=1; i<breaks.size(); ++i){
//        double xsum = 0, Nsum = 0;
//        int count = 0;
//        //cout << "For x < " << x[i] << ":  ";
//        while (pN != state.end() && *pX <= breaks[i]){
//            //cout << *pX << " (" << *pN << "), ";
//            Nsum += *pN;
//            xsum += (*pX)*(*pN);
//            ++count; ++pN; ++pX;
//        }
//        if (Nsum > 0){
//            dens.push_back(Nsum);
//            xmean.push_back(xsum/Nsum);
//        }
//        //cout << endl;
//    }
//    xmean.push_back(xm);

//    vector<float> xbounds(xmean.size()-1); 
//    for (int i=0; i<xbounds.size(); ++i) xbounds[i] = (xmean[i]+xmean[i+1])/2;

//    vector<float> dx(xbounds.size()-1);
//    for (int i=0; i<dx.size(); ++i) dx[i] = xbounds[i+1]-xbounds[i]; 

//    for (int i=0; i<dens.size(); ++i) dens[i] /= dx[i]; 

//    for (auto x : xmean) dens.push_back(x);	
//    return dens;
//}

