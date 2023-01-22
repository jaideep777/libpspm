#ifndef PLANT_PSPM_H_
#define PLANT_PSPM_H_

#include "pspm_environment.h"
#include "plant.h"

#include <solver.h>


class PSPM_Plant : public plant::Plant {
	public:
	
	double t_birth = 0;
	//double input_seed_rain = 1;	

	double viable_seeds;
	double viable_seeds_dt;

	std::vector<std::string> varnames = {"|lma|", "mort", "fec", "vs", "ha", "hm"};
	std::vector<std::string> statevarnames = {"mort", "vs", "ha", "hm"};

	int nrc = 0; // number of evals of compute_vars_phys() - derivative computations actually done by plant
	int ndc = 0; // number of evals of mortality_rate() - derivative computations requested by solver
	int nbc = 0;

	PSPM_Plant(); 
	void set_size(double _x);
	double init_density(double x, void * _env, double input_seed_rain);
	void preCompute(double x, double t, void * _env);
	double establishmentProbability(double t, void * _env);
	double growthRate(double x, double t, void * env);
	double mortalityRate(double x, double t, void * env);
	double birthRate(double x, double t, void * env);
	
	void init_state(double t, void * _env);
	std::vector<double>::iterator set_state(std::vector<double>::iterator &it);
	std::vector<double>::iterator get_state(std::vector<double>::iterator &it);
	std::vector<double>::iterator get_rates(std::vector<double>::iterator &it);
	void print(std::ostream &out = std::cout);

	void save(std::ofstream &fout){}
	void restore(std::ifstream &fin){}

};



#endif 
