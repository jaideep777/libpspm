
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
		dNint[i] = -mod->mortalityRate(xint[i], t);
	}

	double growthGrad = (mod->growthRate(xb+0.001, t) - mod->growthRate(xb, t))/0.001;
	double mortGrad   = (mod->mortalityRate(xb+0.001, t) - mod->mortalityRate(xb, t))/0.001;
	
	double birthFlux = mod->birthRate(x0, t)*N0;
	for (int j=0; j<J; ++j) birthFlux += mod->birthRate(xint[j], t)*Nint[j];

	*dN0  = -mod->mortalityRate(xb, t)*N0 - mortGrad*pi0 + birthFlux;
 	*dpi0 = mod->growthRate(xb, t)*N0 + growthGrad*pi0 - mod->mortalityRate(xb, t)*pi0;
 
}


template <class Model>
void Solver<Model>::addCohort_EBT(){
	auto p_pi0  = state.begin();
	auto p_N0   = advance(state.begin(), J+1); 

	double x0 = xb + *p_pi0/(*p_N0+1e-12); 
	
	state.insert(p_N0,  0); // this inserts 0 BEFORE p_p0. Now both iterators are invalid
	state.insert(state.begin(), 0); // this inserts 0 at the 1st position	

	++J;	// increment J to reflect new system size
}


template <class Model>
void Solver<Model>::removeDeadCohorts_EBT(){
	// since erase() invalidates all iterators beyond erased position, 
	// we flag all elements to be removed in 1 pass, 
	// then in a 2nd pass, iterate backwards removing the flagged cohorts 
	vector<bool> flags(state.size(), false);
	
	auto p_x  = state.begin();
	auto p_N  = advance(state.begin(), J+1); 

	// for dead cohorts, mark N and correponding x for removal  
	++p_x; ++p_N; ++p_x_flag; ++p_N_flag; // skip pi0, N0 from removal
	int Jnew = J;
	for (int i=0; i<J; ++i){
		if (*p_N < 1e-10) *p_x_flag = *p_N_flag = true; 
		----Jnew;
	}
	J = Jnew;

	// remove flagged elements
	auto pend = std::remove_if(state.begin(), state.end(), [&state, &flags](double d){ return flags[&d-&*state.begin()]});
	state.erase(pend, state.end());

}

