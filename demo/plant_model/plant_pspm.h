#ifndef PLANT_PSPM_H_
#define PLANT_PSPM_H_


#include "plant/environment.h"
#include "plant/plant.h"


class LightEnvironment : public plant::Environment{
	//double evalEnv(double x, double t){
	//    env.light_profile.eval(x); // return 1;
	//}

	public:
	LightEnvironment(double openness) : Environment(openness){
	}

	// This function must do any necessary precomputations to facilitate evalEnv()
	// Therefore, this should calculate env for all X when it is a function of X
	// In such a case, the solver's SubdivisionSpline can be ussed
	// Note: The state vector in the solver will not be updated until the RK step is completed. 
	// Hence, explicitly pass the state to this function.
	// ~
	// Also this is the only function that exposes the state vector, so if desired, the state vector 
	// can be saved from here and reused in other rate functions (using createIterators_state())
	// ~
	// TODO: In Solver, add a add_iAttribute() function, that will calculate some individual 
	// level attributes from x, which can be reused if required. E.g., in Plant, we can add leaf_area
	// as an iAttribute. iAttributes can be mapped to integers, say using enums
	// Alternatively, switch to Indiviudual class as a template parameter for solver
	template <class Model>
	void computeEnv(double t, vector<double> &state_vec, Solver<Model, LightEnvironment> * S){
		//            _xm 
		// Calculate / w(z,t)u(z,t)dz
		//        xb`
		auto canopy_openness = [S, t, &state_vec, this](double z){
			double kI = 0.5;

			double leaf_area_above_z = 0;
			
			// Loop over resident species --->
			for (int i=0; i<S->n_species(); ++i){
				plant::Plant * p = &(S->get_species(i)->mod->p);
				auto la_above = [z, p](double x, double t){
					p->set_height(x);	// sets height and leaf-area
					double a = p->area_leaf_above(z, p->vars.height, p->vars.area_leaf);
					return a;	
				};
				leaf_area_above_z += S->integrate_wudx_above(la_above, t, z, state_vec, i);
				//leaf_area_above_z += S->integrate_x(la_above, t, state_vec, i);
			}

			//cout << "LA(" << z << ") = " << exp(-kI*leaf_area_above_z) << "\n";
			return exp(-kI*leaf_area_above_z);
		};	
	
		//cout << S->xb << " " << S->getMaxSize() << endl;	
		time = t;
		//for (int s=0; s<S->n_species(); ++s) S->get_species(s)->u0_save = S->get_u0(t, s);
		light_profile.construct(canopy_openness, 0, S->maxSize(state_vec.begin()));
	}


};


class FixedEnvironment{
	private:
	double openness;

	public:
	FixedEnvironment(double o){
		openness = o;
	}

	double canopy_openness(double z) const {
		return openness;
	}
	
	template <class Model>
	void computeEnv(double t, vector<double> &state_vec, Solver<Model, FixedEnvironment> * S){
	}

	double patch_survival(double t) const{
		return 1;
	}

	
};


class PlantModel{
	public:

	double input_seed_rain = 1;	

	plant::Plant seed; // prototype to be inserted

	int nrc = 0; // number of evals of compute_vars_phys() - derivative computations actually done by plant
	int ndc = 0; // number of evals of mortality_rate() - derivative computations requested by solver
	int nbc = 0;

	// use this to store one-shot rates for each individual and supply through the rate functions
	plant::Plant p;
	
	PlantModel() {
		
	}

	template<class Env>
	double initDensity(double x, Env * env){
		if (x == seed.vars.height){
			p.set_height(x);
			p.compute_vars_phys(*env);
			double u0 = input_seed_rain*p.germination_probability(*env)/growthRate(p.vars.height, 0, env);
			return u0;
		}
		else return 0;
	}


	template<class Env>
	double establishmentProbability(double t, Env * env){
		seed.compute_vars_phys(*env);
		return seed.germination_probability(*env);
	}

	template<class Env>
	double growthRate(double x, double t, Env * env){
		//if (p.vars.height != x){
			p.set_height(x);
			p.compute_vars_phys(*env);
			++nrc;
		//}
		return p.vars.height_dt;
			
	}

	template<class Env>
	double mortalityRate(double x, double t, Env * env){
		//assert(p.vars.height == x);
		if (p.vars.height != x){
			p.set_height(x);
			p.compute_vars_phys(*env);
			++ndc;
		}
		return p.vars.mortality_dt;
	}

	template<class Env>
	double birthRate(double x, double t, Env * env){
		// Need this here because birthRate is not called in order, and only called rarely,
		// after completion of step_to.
		if (p.vars.height != x){
			++nbc;
			p.set_height(x);
			p.compute_vars_phys(*env);
		}
		//assert(p.vars.height == x);
		return p.vars.fecundity_dt;
	}

	
	// optional functions, if extra size-structured variables are desired
	template<class Env>
	vector<double> initStateExtra(double x, double t, Env * env){
		vector<double> sv;
		sv.reserve(4);	
		sv.push_back(-log(seed.germination_probability(*env)/env->patch_survival(t))); // mortality 
		sv.push_back(0); // viable_seeds
		sv.push_back(0); // heartwood area
		sv.push_back(0); // heartwood mass
		return sv;
	}


	template<class Env>
	vector<double>::iterator calcRates_extra(double x, double t, Env * env,
										     vector<double>::iterator istate_ex, vector<double>::iterator irates_ex){
		
		assert(p.vars.height == x);
		
		double p_plant_survival = exp(-(*istate_ex));

		*irates_ex++ = p.vars.mortality_dt;	// mortality
		*irates_ex++ = p.vars.fecundity_dt * env->patch_survival(t) * p_plant_survival; // viable_seeds
		*irates_ex++ = p.vars.area_heartwood_dt; // heartwood area
		*irates_ex++ = p.vars.mass_heartwood_dt; // heartwood mass

		return irates_ex;
	
	}

	//// For output only
	//template<class Env>
	//void setState(Solver<PlantModel, Env> *S){
		
	//    auto iset = S->get_species(0)->get_iterators(S->state);

	//    p.vars.mortality      = *iset.get("mort");
	//    p.vars.fecundity      = *iset.get("fec");
	//    p.vars.area_heartwood = *iset.get("heart_area");
	//    p.vars.mass_heartwood = *iset.get("heart_mass");

	//}


};



#endif 
