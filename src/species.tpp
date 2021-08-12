
template<class Model>
Species<Model>::Species(std::vector<double> breaks){
	J = breaks.size(); // changes with 
	cohorts.resize(J);
	for (int i=0; i<J; ++i) cohorts[i].x = breaks[i];	
}


//template<class Model>
//int Species<Model>::addVar(std::string name, int stride, int offset){
//    varnames.push_back(name);
//    strides.push_back(stride);
//    offsets.push_back(offset);
//}


//template<class Model>
//void Species<Model>::clearVars(){
//    varnames.clear();
//    strides.clear();
//    offsets.clear();
//}


//template<class Model>
//void Species<Model>::set_model(Model *M){
//    mod = M;
//}


void Species_Base::set_inputBirthFlux(double b){
	birth_flux_in = b;
}


//template<class Model>
//double Species<Model>::set_iStateVariables(std::vector<std::string> names){
//    varnames_extra = names;
//}

void Species_Base::set_bfin_is_u0in(bool flag){
	bfin_is_u0in = flag;
}


int Species_Base::xsize(){
	return J;
}


//template<class Model>
//int Species<Model>::size(){
//    return varnames.size()*J;
//}



template<class Model>
double Species<Model>::get_maxSize(){ // TODO ALERT: make sure this sees the latest state
	if (!X.empty()) return *x.rbegin();	// for FMU, get this from X
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


template <class Model>
void Species<Model>::print(){
	std::cout << "~~~~~ Species ~~~~~\n";
	//auto iset = get_iterators(sv);
	//std::cout << "start index = " << start_index <<"\n";
	//std::cout << "Model = " << mod << "\n";
	std::cout << "xsize = " << J << "\n";
	std::cout << "Input birth flux = " << birth_flux_in << "\n";
	//if (!X.empty()){
	//    iset.push_back("_X", X.begin(),1);
	//    iset.push_back("_h", h.begin(),1);
	//    std::cout << "x (" << x.size() << "): "; 
	//    for (auto xx : x) std::cout << xx << " ";
	//    std::cout << "\n";
	//}

	std::cout << "Cohorts:\n";
	for (auto& c : cohorts) std::cout << c.x << " ";
	std::cout << "\n";
	std::cout << "Max size = " << get_maxSize() << "\n";
	//std::cout << "State (" << size() << "):\n";
	//iset.print();
	
	//std::cout << "Rates (" << size() << "):\n";
	//auto irates = get_iterators(rv);
	//irates.print();
	std::cout << "-------\n\n";

}


