#include "pspm_plant.h"

typedef LightEnvironment EnvUsed;


PSPM_Plant::PSPM_Plant() : plant::Plant() {
	
}

void PSPM_Plant::set_size(double _x){
	set_height(_x);
}

double PSPM_Plant::init_density(double x, void * _env, double input_seed_rain){
	//if (x == seed.vars.height){
		//p.set_height(x);
		EnvUsed * env = (EnvUsed*)_env;
		compute_vars_phys(*env);
		double u0;
		if (x == vars.height)
			u0 = input_seed_rain*germination_probability(*env)/vars.height_dt;
		else 
			u0 = 0;
		return u0;
	//}
	//else return 0;
}


void PSPM_Plant::preCompute(double x, double t, void * _env){
	EnvUsed * env = (EnvUsed*)_env;
	compute_vars_phys(*env);
	double p_plant_survival = exp(-vars.mortality);
	//viable_seeds_dt = vars.fecundity_dt; // only for single-plant testrun
	viable_seeds_dt = vars.fecundity_dt * p_plant_survival * env->patch_survival(t) / env->patch_survival(t_birth);
}

double PSPM_Plant::establishmentProbability(double t, void * _env){
	EnvUsed * env = (EnvUsed*)_env;
	return germination_probability(*env);
}

double PSPM_Plant::growthRate(double x, double t, void * env){
	//if (p.vars.height != x){
		//p.set_height(x);
		//compute_vars_phys(*env);
		//++nrc;
	////}
	return vars.height_dt;
		
}

double PSPM_Plant::mortalityRate(double x, double t, void * env){
	//assert(p.vars.height == x);
	//if (p.vars.height != x){
		//p.set_height(x);
		//p.compute_vars_phys(*env);
		//++ndc;
	//}
	return vars.mortality_dt;
}

double PSPM_Plant::birthRate(double x, double t, void * env){
	// Need this here because birthRate is not called in order, and only called rarely,
	// after completion of step_to.
	//if (p.vars.height != x){
		//++nbc;
		//p.set_height(x);
		//p.compute_vars_phys(*env);
	//}
	//assert(p.vars.height == x);
	return vars.fecundity_dt;
}


void PSPM_Plant::init_state(double t, void * _env){
	//set_size(x);	
	EnvUsed * env = (EnvUsed*)_env;
	vars.mortality = -log(germination_probability(*env)); ///env->patch_survival(t));    // mortality
	viable_seeds = 0; // viable seeds
	vars.area_heartwood = 0;   // heartwood area
	vars.mass_heartwood = 0;     // heartwood mass
	t_birth = t;			// set cohort's birth time to current time
	//vars.mortality = 0; // only for single plant testrun
}

vector<double>::iterator PSPM_Plant::set_state(vector<double>::iterator &it){
	vars.mortality      = *it++;
	viable_seeds        = *it++;
	vars.area_heartwood = *it++;
	vars.mass_heartwood = *it++;
	//vars.fecundity = viable_seeds; // only for single plant test run
	return it;
}

vector<double>::iterator PSPM_Plant::get_state(vector<double>::iterator &it){
	*it++ = vars.mortality;
	*it++ = viable_seeds;
	*it++ = vars.area_heartwood;
	*it++ = vars.mass_heartwood;
	return it;
}

vector<double>::iterator PSPM_Plant::get_rates(vector<double>::iterator &it){

	*it++ = vars.mortality_dt;	// mortality
	*it++ = viable_seeds_dt; // viable_seeds
	*it++ = vars.area_heartwood_dt; // heartwood area
	*it++ = vars.mass_heartwood_dt; // heartwood mass
	return it;
}

void PSPM_Plant::print(std::ostream &out){
	out << "|" << lma << "|\t" << vars.mortality << "\t" << vars.fecundity << "\t" << viable_seeds << "\t" << vars.area_heartwood << "\t" << vars.mass_heartwood << "\t";
}

