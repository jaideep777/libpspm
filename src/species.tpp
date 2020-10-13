


template<class Model>
int Species<Model>::addVar(std::string name, int stride, int offset){
	varnames.push_back(name);
	strides.push_back(stride);
	offsets.push_back(offset);
}


template<class Model>
void Species<Model>::clearVars(){
	varnames.clear();
	strides.clear();
	offsets.clear();
}


template<class Model>
void Species<Model>::set_model(Model *M){
	mod = M;
}


template<class Model>
void Species<Model>::set_inputBirthFlux(double b){
	birth_flux_in = b;
}


template<class Model>
double Species<Model>::set_iStateVariables(std::vector<std::string> names){
	varnames_extra = names;
}


template<class Model>
int Species<Model>::xsize(){
	return J;
}


template<class Model>
double Species<Model>::get_maxSize(std::vector<double>::iterator state_begin){
	if (!X.empty()) return *x.rbegin();
	else {
		return *next(state_begin, start_index + J-1);
	}
}


template<class Model>
IteratorSet<std::vector<double>::iterator> Species<Model>::get_iterators(std::vector<double> &v){
	IteratorSet<std::vector<double>::iterator> iset(v.begin()+start_index, varnames, J, offsets, strides);
	if (!X.empty()) iset.push_back("X", X.begin(), 1);
	return iset;
}


