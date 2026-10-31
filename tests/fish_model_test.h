#include <iostream>
#include <cmath>
#include <environment_base.h>
#include <individual_base.h>
using namespace std;



class SeaEnv : public EnvironmentBase{
	public:
	void computeEnv(double t, Solver * S, vector<double>::iterator s, vector<double>::iterator dsdt){
	}

};

const int DIM = 1, CAT = 2;

class Fish : public IndividualBase<DIM, CAT>{
	public:
	double init_density(void * _env, double bf){
		return 1.0/x[0];
	}

	void set_size(const array<double,DIM>& x){
	}

	array<double,DIM> growthRate(double t, void * env){
		return {1};
	}
	double mortalityRate(double t, void * env){
		return 0.5;
	}
	double birthRate(double t, void * env){
		return 0.5;
	}

	size_t init_category(double t, void * _env){
		return 0;
	}

	void set_category(size_t _q){
	}

	std::vector<std::vector<double>> get_transition_matrix(double t, double dt, void * _env){
		return {{0.8*dt, 0.2*dt}, {0.0*dt, 1.0*dt}};
	}

};

