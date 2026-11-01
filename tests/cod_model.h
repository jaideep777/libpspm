#include <iostream>
#include <cmath>
#include <environment_base.h>
#include <individual_base.h>
using namespace std;

class SeaEnv;

const int DIM = 1, CAT = 2;


struct FishParams {
	// biology
	double theta = 8.10e-6; // kg/cm 
	double zeta  = 3.01; 
	
	double amax = 30;

	std::vector<double> ma = {0,  0.0000,  0.0000,  0.0002,  0.0029,  0.0565,  0.3240,  0.6361,  0.8451,  0.9473,  0.9838,  0.9944,  0.9964,  1.0000,  1.0000,  1.0000, 1.0000,  1.0000,  
		                          1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000};   // maturation probability with age
	std::vector<double> mai = {0, 1.48936, 0.61272, 0.36018, 0.26108, 0.22300, 0.20447, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000,
								  0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000}; // mortality rate of immature indivuduals
	std::vector<double> mam = {0, 1.48936, 0.61272, 0.36018, 0.26108, 0.22300, 0.20447, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000,
								  0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000}; // mortality rate of mature indivuduals

	double l8 = 196;
	double kappa = 0.0555;
	double a0 = -0.937;

	double r0 = 21.77072;		// recruitment rate per kg SSB 
	double Bhalf = 365426284;	// Half saturation constant of recruitment

	// debug
	int flag = 0;
};


class Fish{

	public:
	FishParams par;
	
	// state variables
	double age = 1;		// age in years
	double length;  // length in cm

	// physiological variables
	double weight;	// weight in kg

	public:
	void set_age(int _a){
		age = _a;

		// in an explicity age-length model, this will go.
		length = par.l8*(1-exp(-par.kappa*(age-par.a0)));
		set_length(length);	
	}

	void set_length(double s){
		length = s;
		weight = par.theta*pow(length, par.zeta);  // Eq. 2
	}
	
	double maturation_rate(){
		int a  = int(age);
		return (a > par.amax)? 1e20 : (par.ma[a+1]-par.ma[a])/(1 - par.ma[a]);
	}
	
	double naturalMortalityRate(bool isMature){
		int a  = int(age);
		double rate;
		if (isMature) rate = (a > par.amax)? 1e20 : par.mam[a]; // FIXME: use inf
		else          rate = (a > par.amax)? 1e20 : par.mai[a]; // FIXME: use inf
		return rate;
	}

};



class PSPMFish : public IndividualBase<DIM, CAT>, public Fish{
	public:
	double init_density(void * _env, double bf){
		return 1.0e9/x[0];
	}

	// "size" is age in this model
	void set_size(const array<double,DIM>& x){
		Fish::set_age(x[0]);
	}

	// age-growth rate is always 1
	array<double,DIM> growthRate(double t, void * env){
		return {1};
	}

	double mortalityRate(double t, void * env){
		return Fish::naturalMortalityRate(q);
	}

	double birthRate(double t, void * _env);

	size_t init_category(double t, void * _env){
		return 0;
	}

	void set_category(size_t _q){
	}

	std::vector<std::vector<double>> get_transition_matrix(double t, void * _env){
		double p = Fish::maturation_rate();
		p = (p > 1)? 1 : p; // cap at 1
		return {{1-p, p}, {0.0, 1.0}};
		// return {{0.8, 0.2}, {0.0, 1.0}};
	}

};

class SeaEnv : public EnvironmentBase{
	public:
	double ssb = 0, tsb = 0;

	void computeEnv(double t, Solver * S, vector<double>::iterator s, vector<double>::iterator dsdt){
		auto lambda_ssb = [S](int i, double t){
			auto f = static_cast<Species<PSPMFish>*>(S->species_vec[0])->getCohort(i);
			PSPMFish * fish = static_cast<PSPMFish*>(&f);
			return fish->q * fish->weight;
		};
		ssb = S->state_integral(lambda_ssb, t, 0);

		auto lambda_tsb = [S](int i, double t){
			auto f = static_cast<Species<PSPMFish>*>(S->species_vec[0])->getCohort(i);
			PSPMFish * fish = static_cast<PSPMFish*>(&f);
			return fish->weight;
		};
		tsb = S->state_integral(lambda_tsb, t, 0);

	}

};


double PSPMFish::birthRate(double t, void * _env){
	auto env = static_cast<SeaEnv*>(_env);
	return q*par.r0*weight / (1 + env->ssb/par.Bhalf);
	// return 0;
}
