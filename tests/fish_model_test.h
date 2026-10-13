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

const int DIM = 1;

class Fish : public IndividualBase<DIM>{
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
		return 0.1;
	}

};

