#ifndef  PSPM_PSPM_COHORT_H_
#define  PSPM_PSPM_COHORT_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include "io_utils.h"

template<class Ind>
class Cohort : public Ind {
	public:
	static int np, ng, nm, nf;   // number of evaluations of demographic functions

	double x = -999;
	std::vector<double> xn = { -999 }; // is this initialisation needed?
	double u = -999;
	int id = 0;	

	double birth_time = 0;
	bool remove = false;
	
	bool need_precompute = true;

	// Default constructor simply calls Individual's default ctor
	// NOTE: This means user's Ind class will need ctor with no arguments
	Cohort() : Ind() {
	}

	// Construct a cohort from Individual using copy constructor of Individual
	Cohort(const Ind& _ind) : Ind(_ind){
	}

	void print_xu(std::ostream &out = std::cout){
		out << std::setw(6)  << std::setprecision(4) << birth_time 
		    << std::setw(12) << std::setprecision(4) << x 
			<< std::setw(12) << std::setprecision(4) << u; 
	}

	void print(std::ostream &out = std::cout){
		print_xu(out);
		Ind::print(out);
	}

	// This function is temporary until fully transitioned all solvers to multi-state format
	void set_size(double _x){
		x = _x;
		std::vector<double> _xn;
    	_xn.push_back(x);
		set_size(_xn);
	}

	void set_size(std::vector <double> _xn){
		xn = _xn;
		need_precompute = true;
		Ind::set_size(_xn);
	}
	
	
	//  These are defined here so that precompute trigger can be 
	//  checked before calling user-defined function 
	void preCompute(double x, double t, void * _env){
		++np;
		//std::cout << "cohort precompute: "; print(); std::cout << "\n";
		Ind::preCompute(x,t,_env);	
		need_precompute = false;   // Once precompute is called, no need to further precompute until necessary
	}

	void preCompute(std::vector <double> xn, double t, void * _env){
		++np;
		//std::cout << "cohort precompute: "; print(); std::cout << "\n";
		Ind::preCompute(xn,t,_env);	
		need_precompute = false;   // Once precompute is called, no need to further precompute until necessary
	}
	
	double growthRate(double x, double t, void * _env){
		++ng;
		if (need_precompute) preCompute(x,t,_env);
		//std::cout << "cohort growthRate(): "; print(); std::cout << "\n";
		return Ind::growthRate(x,t,_env);	
	}

	std::vector<double> growthRate(std::vector<double> xn, double t, void * _env){
		++ng;
		if (need_precompute) preCompute(xn,t,_env);
		//std::cout << "cohort growthRate(): "; print(); std::cout << "\n";
		return Ind::growthRate(xn,t,_env);	
	}
	
	double mortalityRate(double x, double t, void * _env){
		++nm;
		if (need_precompute) preCompute(x,t,_env);
		//std::cout << "cohort mortRate(): "; print(); std::cout << "\n";
		return Ind::mortalityRate(x,t,_env);	
	}
	
	double mortalityRate(std::vector <double> xn, double t, void * _env){
		++nm;
		if (need_precompute) preCompute(xn,t,_env);
		//std::cout << "cohort mortRate(): "; print(); std::cout << "\n";
		return Ind::mortalityRate(xn,t,_env);	
	}

	double birthRate(double x, double t, void * _env){
		++nf;
		if (need_precompute) preCompute(x,t,_env);
		//std::cout << x << " cohort birthRate: "; print(); std::cout << "\n";
		return Ind::birthRate(x,t,_env);	
	}
	
	double birthRate(std::vector<double> xn, double t, void * _env){
		++nf;
		if (need_precompute) preCompute(xn,t,_env);
		//std::cout << x << " cohort birthRate: "; print(); std::cout << "\n";
		return Ind::birthRate(xn,t,_env);	
	}

	void save(std::ofstream &fout, int n_extra_vars){
		// Save/restore individual first (for metadata)
		Ind::save(fout);

		// Then save cohort (for cohort state). This way, Individual metadata will be available when set_size() and set_state() are called in restore()
		fout << "Cohort<Ind>::v1" << "   ";
		fout << std::make_tuple(
				  id
				, birth_time
				, remove
				, need_precompute); // we actually need not save need_precompute, because set_size() will always set it to 1 during restore
		fout << x << ' ' << u << ' '; // TODO Edit to to make it multistate 

		std::vector<double> ex_state(n_extra_vars);
		auto it = ex_state.begin();
		Ind::get_state(it);
		fout << ex_state;
	}

	void restore(std::ifstream &fin, int n_extra_vars){
		Ind::restore(fin);

		std::string s; fin >> s; // discard version number
		fin >> id
		    >> birth_time
		    >> remove
		    >> need_precompute;

		fin >> x >> u; // TODO Edit to to make it multistate 
		set_size(x);

		std::vector<double> ex_state(n_extra_vars);
		fin >> ex_state;
		auto it = ex_state.begin();
		Ind::set_state(it);
	}
	
	// FIXME other env dependent rates should also check for precompute
};

template<class Model>
int Cohort<Model>::np = 0;

template<class Model>
int Cohort<Model>::ng = 0;

template<class Model>
int Cohort<Model>::nm = 0;

template<class Model>
int Cohort<Model>::nf = 0;

#endif

