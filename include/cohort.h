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

	// std::array<double,dim> x; // Moved to IndividualBase
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
		out << std::setw(6)  << std::setprecision(4) << birth_time;
		for (auto xx : Ind::x) out << std::setw(12) << std::setprecision(4) << xx; 
		out << std::setw(12) << std::setprecision(4) << u; 
	}

	void print(std::ostream &out = std::cout){
		print_xu(out);
		Ind::print(out);
	}

	template<size_t dim>
	void set_size(std::array<double,dim> _x){
		Ind::x = _x;
		need_precompute = true;
		Ind::set_size(_x);
	}
	
	
	//  These are defined here so that precompute trigger can be 
	//  checked before calling user-defined function 
	template<size_t dim>
	void preCompute(std::array<double,dim> x, double t, void * _env){
		++np;
		//std::cout << "cohort precompute: "; print(); std::cout << "\n";
		Ind::preCompute(x,t,_env);	
		need_precompute = false;   // Once precompute is called, no need to further precompute until necessary
	}
	

	template<size_t dim>
	std::array<double,dim> growthRate(std::array<double,dim> x, double t, void * _env){
		++ng;
		if (need_precompute) preCompute(x,t,_env);
		//std::cout << "cohort growthRate(): "; print(); std::cout << "\n";
		return Ind::growthRate(x,t,_env);	
	}
	

	template<size_t dim>
	double mortalityRate(std::array<double,dim> x, double t, void * _env){
		++nm;
		if (need_precompute) preCompute(x,t,_env);
		//std::cout << "cohort mortRate(): "; print(); std::cout << "\n";
		return Ind::mortalityRate(x,t,_env);	
	}
	

	template<size_t dim>
	double birthRate(std::array<double,dim> x, double t, void * _env){
		++nf;
		if (need_precompute) preCompute(x,t,_env);
		//std::cout << x << " cohort birthRate: "; print(); std::cout << "\n";
		return Ind::birthRate(x,t,_env);	
	}
	

	void save(std::ostream &fout, int n_extra_vars){
		// Save/restore individual first (for metadata)
		Ind::save(fout);

		// Then save cohort (for cohort state). This way, Individual metadata will be available when set_size() and set_state() are called in restore()
		fout << "Cohort<Ind>::v1" << "   ";
		fout << std::make_tuple(
				  id
				, birth_time
				, remove
				, need_precompute); // we actually need not save need_precompute, because set_size() will always set it to 1 during restore
		fout << Ind::x;
		fout << u << ' ';

		std::vector<double> ex_state(n_extra_vars);
		auto it = ex_state.begin();
		Ind::get_state(it);
		fout << ex_state;
	}


	void restore(std::istream &fin, int n_extra_vars){
		Ind::restore(fin);

		std::string s; fin >> s; // discard version number
		fin >> id
		    >> birth_time
		    >> remove
		    >> need_precompute;

		fin >> Ind::x;
		fin >> u;
		set_size(Ind::x);

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

