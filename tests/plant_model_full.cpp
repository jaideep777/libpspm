#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
#include "solver.h"

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
			}

			//cout << "la = " << leaf_area_above << "\n";
			return exp(-kI*leaf_area_above_z);
		};	
	
		//cout << S->xb << " " << S->getMaxSize() << endl;	
		time = t;
		light_profile.construct(canopy_openness, 0, S->maxSize(state_vec.begin()));
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

	double initDensity(double x, LightEnvironment * env){
		if (x == seed.vars.height){
			p.set_height(x);
			p.compute_vars_phys(*env);
			double u0 = input_seed_rain*p.germination_probability(*env)/growthRate(p.vars.height, 0, env);
			return u0;
		}
		else return 0;
	}



	double establishmentProbability(double t, LightEnvironment * env){
		seed.compute_vars_phys(*env);
		return seed.germination_probability(*env);
	}

	double growthRate(double x, double t, LightEnvironment * env){
		//if (p.vars.height != x){
			p.set_height(x);
			p.compute_vars_phys(*env);
			++nrc;
		//}
		return p.vars.height_dt;
			
	}

	double mortalityRate(double x, double t, LightEnvironment * env){
		assert(p.vars.height == x);
		++ndc;
		return p.vars.mortality_dt;
	}

	double birthRate(double x, double t, LightEnvironment * env){
		// Need this here because birthRate is not called in order, and only called rarely,
		// after completion of step_to.
		if (p.vars.height != x){
			++nbc;
			p.set_height(x);
			p.compute_vars_phys(*env);
			++nrc;
		}
		//assert(p.vars.height == x);
		return p.vars.fecundity_dt;
	}

	
	// optional functions, if extra size-structured variables are desired
	vector<double> initStateExtra(double x, double t, LightEnvironment * env){
		vector<double> sv;
		sv.reserve(4);	
		sv.push_back(-log(seed.germination_probability(*env)/env->patch_survival(t))); // mortality 
		sv.push_back(0); // viable_seeds
		sv.push_back(0); // heartwood area
		sv.push_back(0); // heartwood mass
		return sv;
	}


	vector<double>::iterator calcRates_extra(double x, double t, LightEnvironment * env,
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
	//void setState(Solver<PlantModel, LightEnvironment> *S){
		
	//    auto iset = S->get_species(0)->get_iterators(S->state);

	//    p.vars.mortality      = *iset.get("mort");
	//    p.vars.fecundity      = *iset.get("fec");
	//    p.vars.area_heartwood = *iset.get("heart_area");
	//    p.vars.mass_heartwood = *iset.get("heart_mass");

	//}


};




vector<double> generateDefaultCohortSchedule(double max_time){

	vector<double> tvec;

	const double multiplier=0.2, min_step_size=1e-5, max_step_size=2.0;
	
	assert(min_step_size > 0 && "The minimum step size must be greater than zero");
	
	double dt = 0.0, time = 0.0;
	tvec.push_back(time);
	while (time <= max_time) {
		dt = exp2(floor(log2(time * multiplier)));
		time += min(max(dt, min_step_size), max_step_size);
		tvec.push_back(time);
	}

	// Drop the last time; that's not going to be needed:
	if (tvec.size() >=1) 	// JAI: added to avoid overflow warning
		tvec.resize(tvec.size() - 1);

	return tvec;
}


template<class Model, class Environment>
class SolverIO{
	public:
	int nspecies;
	Solver<Model, Environment> * S;

	vector <vector<ofstream>> streams;

	void openStreams(){
		for (int s=0; s < S->n_species(); ++s){
			auto spp = S->get_species(s);
			vector<string> varnames = spp->get_varnames();
			vector<ofstream> spp_streams;
			for (string name : varnames){
				stringstream sout;
				sout << "species_" << s << "_" << name + ".txt";
				cout << sout.str() << endl;
				ofstream fout(sout.str().c_str());
				spp_streams.push_back(std::move(fout));
			}
			streams.push_back(std::move(spp_streams));
		}
	}

	void closeStreams(){
		for (int s=0; s<streams.size(); ++s){
			for (int j=0; j<streams[s].size(); ++j){
				streams[s][j].close();
			}
		}
	}

	void writeState(){
		for (int s=0; s < S->n_species(); ++s){
			auto spp = S->get_species(s);
			auto iset = spp->get_iterators(S->state);
			auto& ivec = iset.get();

			for (int i=0; i<streams[s].size(); ++i) streams[s][i] << S->current_time << "\t";

			for (iset.rbegin(); !iset.rend(); --iset){
				for (int i=0; i<ivec.size(); ++i){
					streams[s][i] << *ivec[i] << "\t"; 
				}
			}
			
			for (int i=0; i<streams[s].size(); ++i) streams[s][i] << "\n";
		}
	}
};


int main(){
	
	//initPlantParameters(plant::par);
	
	LightEnvironment env(1);	
	env.light_profile.print();	
	
	//plant::Environment env(1);
	plant::Plant p1;

	plant::Plant p;
	p.lma = 0.2625;
	p.initParameters();
	p.vars.height = p.par.height_0; //0.3257146; //0.3920458; //0.3441948;
	p.vars.area_leaf = p.par.area_leaf_0; 


	plant::Plant p3;
	p3.lma = 0.4625;
	p3.initParameters();
	p3.vars.height = p3.par.height_0; //0.3257146; //0.3920458; //0.3441948;
	p3.vars.area_leaf = p3.par.area_leaf_0; 

	
	cout << p1 << endl;
	cout << p << endl;

	//exit(1);

    Solver<PlantModel, LightEnvironment> S(SOLVER_CM);
    S.use_log_densities = true;
	S.control.ode_eps = 1e-4;
	S.setEnvironment(&env);
	//    S.createSizeStructuredVariables({"mort", "fec", "heart_area", "heart_mass"});
    
	PlantModel M1;
	M1.p = M1.seed = p1;
    cout << "HT1 === " << M1.p.vars.height << endl;
	

    PlantModel M;
    M.p = M.seed = p;
    cout << "HT === " << M.p.vars.height << endl;


    PlantModel M3;
    M3.p = M3.seed = p3;
    cout << "HT === " << M3.p.vars.height << endl;

	S.addSpecies(vector<double>(1, M1.p.vars.height), &M1, {"mort", "fec", "heart", "sap"}, M1.input_seed_rain);
	S.addSpecies(vector<double>(1, M.p.vars.height), &M, {"mort", "fec", "heart", "sap"}, M.input_seed_rain);
	S.addSpecies(vector<double>(1, M3.p.vars.height), &M3, {"mort", "fec", "heart", "sap"}, M3.input_seed_rain);
	
	S.resetState();
    S.initialize();

    S.print();

	vector <double> times = generateDefaultCohortSchedule(105.32);
	for (auto t : times) cout << t << " "; cout << endl;

	
	SolverIO<PlantModel, LightEnvironment> sio;
	sio.S = &S;
	sio.openStreams();

	ofstream fli("light_profile_ind_plant.txt");
	
	vector <vector<double>> seeds_out(S.n_species());

	for (size_t i=0; i < times.size(); ++i){

		S.step_to(times[i]);		
		
		vector<double> seeds = S.newborns_out();
		for (int s=0; s< S.n_species(); ++s){
			double S_D = 0.25;
			seeds_out[s].push_back(seeds[s] * S_D * env.patch_age_density(times[i]));
		}

		cout << times[i] << " " << S.get_species(0)->xsize() << " " << env.light_profile.npoints << " | " << M.nrc << " " << M.ndc << " | " << M.nbc <<"\n";

		vector<double> xl = seq(0, 20, 200);
		for (auto h : xl) fli << env.canopy_openness(h) << "\t";
		fli << endl;

		sio.writeState();

	}
	
	fli.close();
	sio.closeStreams();
	cout << "derivative computations requested/done: " << M.nrc << " " << M.ndc << endl;

	for (int s=0; s< S.n_species(); ++s){
		auto spp = S.get_species(s);
		auto iset = spp->get_iterators(S.state);
		auto& itf = iset.get("fec");
		vector <double> fec_vec;
		fec_vec.reserve(spp->xsize());
		iset.rbegin();
		for (int i=0; !iset.rend(); --iset, ++i){
			double patch_age_density = env.patch_age_density(times[i]);
			double S_D = 0.25;
			double output_seeds = spp->mod->input_seed_rain * S_D * patch_age_density * (*itf);
			//cout << times[i] << " " << M.input_seed_rain << " " << S_D << " " << patch_age_density << " " << (*itf) << " | " << output_seeds << endl;
			fec_vec.push_back(output_seeds);
		}
		cout << "Seed rain for Species " << s << " = " << pn::integrate_trapezium(times, fec_vec) << endl;
		cout << "Seed rain for Species " << s << " (new method) = " << pn::integrate_trapezium(times, seeds_out[s]) << endl;

	}

}

