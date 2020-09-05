#include <algorithm>
#include <cassert>

template <class Model>
void Solver<Model>::calcRates_EBT(double t, vector<double>&S, vector<double> &dSdt){
	double   pi0  =  S[0];
	double * xint = &S[1];
	double   N0   =  S[J+1];
	double * Nint = &S[J+2];

	double * dpi0  = &dSdt[0];
	double * dxint = &dSdt[1];
	double * dN0   = &dSdt[J+1];
	double * dNint = &dSdt[J+2];

	double x0 = xb + pi0/(N0+1e-12); 

	for (size_t i=0; i<J; ++i){
		dxint[i] =  mod->growthRate(xint[i], t);
		dNint[i] = -mod->mortalityRate(xint[i], t)*Nint[i];
	}

	double grad_dx = 0.001;
	double growthGrad = (mod->growthRate(xb+0.001, t) - mod->growthRate(xb, t))/0.001;
	double mortGrad   = (mod->mortalityRate(xb+0.001, t) - mod->mortalityRate(xb, t))/0.001;
	
	double birthFlux = integrate_x([this](double z, double t){return mod->birthRate(z,t);}, t, S, 1);
	
	double B = mod->birthRate(x0, t)*N0;
	for (int j=0; j<J; ++j) B += mod->birthRate(xint[j], t)*Nint[j];
	//cout << "t/b = " << t << " " << mod->evalEnv(0,t) << " " << birthFlux << " ";
	//cout << bf << endl;
	assert(B == birthFlux);

	*dN0  = -mod->mortalityRate(xb, t)*N0 - mortGrad*pi0 + birthFlux;
 	*dpi0 = mod->growthRate(xb, t)*N0 + growthGrad*pi0 - mod->mortalityRate(xb, t)*pi0;
 
}


template <class Model>
void Solver<Model>::addCohort_EBT(){
	auto p_pi0  = state.begin();
	auto p_N0   = state.begin(); advance(p_N0, J+1); 

	double x0 = xb + *p_pi0/(*p_N0+1e-12); 
	
	*p_pi0 = x0; // pi0 becomes x0 when boundary cohort is internalized

	state.insert(p_N0,  0); // insert new N0 (= 0) BEFORE p_p0. (Now both iterators are invalid)
	state.insert(state.begin(), 0); // this inserts new pi0 (= 0) at the 1st position	

	++J;	// increment J to reflect new system size
}


// this function was tested here: cpp.sh/67t3i
template <class Model>
void Solver<Model>::removeDeadCohorts_EBT(){
	// since erase() invalidates all iterators beyond erased position, 
	// we flag all elements to be removed in 1 pass, 
	// then in a 2nd pass, iterate backwards removing the flagged cohorts 
	// One could have simply set the values of flagged x and N to some random value (e.g. 1e-20) and use it in remove_if(); but that is not fully failsafe
	// using iterators so that this code works even if state is a list
	vector<bool> flags(state.size(), false);
	
	// iterators to beginning of x and beginning of N in state array
    auto p_N = state.begin(); advance(p_N, J+1);
	auto p_x_flag = flags.begin();
	auto p_N_flag = flags.begin(); advance(p_N_flag, J+1);

	// for dead cohorts, mark N and correponding x for removal  
	++p_N; ++p_x_flag; ++p_N_flag; // skip pi0, N0 from removal
	int Jnew = J;
	for (int i=0; i<J; ++i){
		if (*p_N < 1e-10){
			*p_x_flag = *p_N_flag = true; 
			--Jnew;
		}
		++p_N; ++p_x_flag; ++p_N_flag; 
	}
	J = Jnew;

	// remove flagged elements
	auto pred = [this, &flags](const double& s) -> bool {
		return flags[&s - (const double*)&this->state[0]]; // TODO: make conatiner-type safe
	};

	auto pend = std::remove_if(state.begin(), state.end(), pred);
	state.erase(pend, state.end());

}



// FIXME: Need to be fixed. Cross check with R
template <class Model>
vector<double> Solver<Model>::cohortsToDensity_EBT(vector <double> &breaks){

	auto pX = state.begin();
	auto pN = state.begin(); advance(pN, J+1);

	*pX = xb + *pX/(*pN+1e-12); 
				
	vector<double> dens;
	vector<double> xmean(1, xb);

	for (int i=1; i<breaks.size(); ++i){
		double xsum = 0, Nsum = 0;
		int count = 0;
		//cout << "For x < " << x[i] << ":  ";
		while (pN != state.end() && *pX <= breaks[i]){
			//cout << *pX << " (" << *pN << "), ";
			Nsum += *pN;
			xsum += (*pX)*(*pN);
			++count; ++pN; ++pX;
		}
		if (Nsum > 0){
			dens.push_back(Nsum);
			xmean.push_back(xsum/Nsum);
		}
		//cout << endl;
	}
	xmean.push_back(xm);

	vector<float> xbounds(xmean.size()-1); 
	for (int i=0; i<xbounds.size(); ++i) xbounds[i] = (xmean[i]+xmean[i+1])/2;

	vector<float> dx(xbounds.size()-1);
	for (int i=0; i<dx.size(); ++i) dx[i] = xbounds[i+1]-xbounds[i]; 

	for (int i=0; i<dens.size(); ++i) dens[i] /= dx[i]; 

	for (auto x : xmean) dens.push_back(x);	
	return dens;
}

