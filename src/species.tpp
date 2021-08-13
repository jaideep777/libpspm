
template<class Model>
void Species<Model>::resize(int _J){
	J = _J;
	cohorts.resize(J);
}


template<class Model>
double Species<Model>::get_maxSize(){ // TODO ALERT: make sure this sees the latest state
	if (!X.empty()) return *x.rbegin();	// for FMU, get this from X
	else if (cohorts.empty()) return 0;
	else {								// else get from state vector
		return cohorts[J-1].x;
	}
}


//template<class Model>
//std::vector<std::string> Species<Model>::get_varnames(){
//    return varnames;
//}


//template<class Model>
//IteratorSet<std::vector<double>::iterator> Species<Model>::get_iterators(std::vector<double> &v){
//    IteratorSet<std::vector<double>::iterator> iset(v.begin()+start_index, varnames, J, offsets, strides);
//    if (!X.empty()) iset.push_back("X", X.begin(), 1);
//    return iset;
//}



template<class Model>
Species<Model>::Species(std::vector<double> breaks){
	J = breaks.size(); // changes with 
	cohorts.resize(J);
	for (int i=0; i<J; ++i) cohorts[i].x = breaks[i];	
}


template <class Model>
void Species<Model>::print(){
	std::cout << "~~~~~ Species ~~~~~\n";
	//auto iset = get_iterators(sv);
	std::cout << "start index = " << start_index <<"\n";
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

	std::cout << "Cohorts:\n";
	std::cout << "id" << "\tx" << "\t" << "u" << "\t";
	if (!cohorts.empty()){
		for (auto s : cohorts[1].varnames) std::cout << s << "\t";
	}
	std::cout << "\n";
	for (auto& c : cohorts) {
		c.print_xu(); //std::cout << c.x << "\t" << c.u << "\t";
		c.print();
		std::cout << "\n";
	}
	std::cout << "\n";
	std::cout << "Max size = " << get_maxSize() << "\n";
	//std::cout << "State (" << size() << "):\n";
	//iset.print();
	
	//std::cout << "Rates (" << size() << "):\n";
	//auto irates = get_iterators(rv);
	//irates.print();
	std::cout << "-------\n\n";

}


template <class Model>
double Species<Model>::getX(int i){
	return cohorts[i].x;
}


template <class Model>
double Species<Model>::getU(int i){
	return cohorts[i].u;
}


template <class Model>
void Species<Model>::set_xb(double _xb){
	xb = _xb;
	boundaryCohort.x = _xb;
}


template <class Model>
void Species<Model>::setX(int i, double _x){
	cohorts[i].x = _x;
}


template <class Model>
void Species<Model>::setU(int i, double _u){
	cohorts[i].u = _u;
}


template <class Model>
void Species<Model>::init_ExtraState(std::vector<double>::iterator &it){
	for (auto& c : cohorts) c.init_state(c.x, it);
}


template <class Model>
double Species<Model>::init_density(int i, double _x){
	return cohorts[i].init_density(_x);
}


template <class Model>
void Species<Model>::copyExtraStateToCohorts(std::vector<double>::iterator &it){
	for (auto& c : cohorts) c.set_state(it);
}


template <class Model>
double Species<Model>::get_u0(double t, void * env){
	
	if (bfin_is_u0in){
		boundaryCohort.u = birth_flux_in;
	}
	else {	
		double g = boundaryCohort.growthRate(boundaryCohort.x, t, env); 
		boundaryCohort.u = (g>0)? birth_flux_in * boundaryCohort.establishmentProbability(t, env)/g  :  0; 
	}
	return boundaryCohort.u;
}

