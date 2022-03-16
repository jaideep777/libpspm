#ifndef  PSPM_PSPM_COHORT_H_
#define  PSPM_PSPM_COHORT_H_

#include <iostream>

template<class Ind>
class Cohort : public Ind {
	public:
	static int np, ng, nm, nf;   // number of evaluations of demographic functions

	double x = -999;
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
		out << birth_time << "\t" << x << "\t" << u << "\t"; 
	}

	void print(std::ostream &out = std::cout){
		print_xu(out);
		Ind::print(out);
	}

	void set_size(double _x){
		x = _x;
		need_precompute = true; // when size is updated, next rate calc will need precompute
		Ind::set_size(x);
	}
	
	
	//  These are defined here so that precompute trigger can be 
	//  checked before calling user-defined function 
	void preCompute(double x, double t, void * _env){
		++np;
		//std::cout << "cohort precompute: "; print(); std::cout << "\n";
		Ind::preCompute(x,t,_env);	
		need_precompute = false;   // Once precompute is called, no need to further precompute until necessary
	}
	
	double growthRate(double x, double t, void * _env){
		++ng;
		if (need_precompute) preCompute(x,t,_env);
		//std::cout << "cohort growthRate(): "; print(); std::cout << "\n";
		return Ind::growthRate(x,t,_env);	
	}
	
	double mortalityRate(double x, double t, void * _env){
		++nm;
		if (need_precompute) preCompute(x,t,_env);
		//std::cout << "cohort mortRate(): "; print(); std::cout << "\n";
		return Ind::mortalityRate(x,t,_env);	
	}
	
	double birthRate(double x, double t, void * _env){
		++nf;
		if (need_precompute) preCompute(x,t,_env);
		//std::cout << x << " cohort birthRate: "; print(); std::cout << "\n";
		return Ind::birthRate(x,t,_env);	
	}
	
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

