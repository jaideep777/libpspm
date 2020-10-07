#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
#include "solver.h"
#include "plant/plant.h"

class PlantModel{
	public:

	double input_seed_rain = 200;	
	plant::Environment env;

	plant::Plant seed; // prototype to be inserted

	int nrc = 0; // number of evals of compute_vars_phys() - derivative computations actually done by plant
	int ndc = 0; // number of evals of mortality_rate() - derivative computations requested by solver

	// use this to store one-shot rates for each individual and supply through the rate functions
	plant::Plant p;
	
	PlantModel() : env(1) {
		
	}

	double initDensity(double x){
		if (x == seed.vars.height){
			seed.set_height(x);
			seed.compute_vars_phys(env);
			double u0 = input_seed_rain*p.germination_probability(env)/growthRate(p.vars.height, 0);
			return u0;
		}
		else return 0;
	}


	//double evalEnv(double x, double t){
	//    env.light_profile.eval(x); // return 1;
	//}
	
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
	void computeEnv(double t, vector<double> &state_vec, Solver<PlantModel> * S){
		//            _xm 
		// Calculate / w(z,t)u(z,t)dz
		//        xb`
		auto canopy_openness = [S, t, &state_vec, this](double z){
			double kI = 0.5;

			auto la_above = [z, this](double x, double t){
				p.set_height(x);
				double a = p.area_leaf_above(z, p.vars.height, p.vars.area_leaf);
				return a;	
			};
			double leaf_area_above = S->integrate_wudx_above(la_above, t, z, state_vec);
			//cout << "la = " << leaf_area_above << "\n";
			return exp(-kI*leaf_area_above);
		};	
	
		//cout << S->xb << " " << S->getMaxSize() << endl;	
		env.time = t;
		env.light_profile.construct(canopy_openness, 0, S->getMaxSize(state_vec.begin()));
	}


	double establishmentProbability(double t){
		seed.compute_vars_phys(env);
		return seed.germination_probability(env);
	}

	double growthRate(double x, double t){
		//if (p.vars.height != x){
			p.set_height(x);
			p.compute_vars_phys(env);
			++nrc;
		//}
		return p.vars.height_dt;
			
	}

	double mortalityRate(double x, double t){
		assert(p.vars.height == x);
		++ndc;
		return p.vars.mortality_dt;
	}

	double birthRate(double x, double t){
		assert(p.vars.height == x);
		return p.vars.fecundity_dt;
	}

	
	// optional functions, if extra size-structured variables are desired
	vector<double> initStateExtra(double x, double t){
		vector<double> sv;
		sv.reserve(4);	
		sv.push_back(-log(seed.germination_probability(env)/env.patch_survival(t))); // mortality 
		sv.push_back(0); // viable_seeds
		sv.push_back(0); // heartwood area
		sv.push_back(0); // heartwood mass
		return sv;
	}


	vector<double>::iterator calcRates_extra(double t, double x, vector<double>::iterator istate_ex, vector<double>::iterator irates_ex){
		
		assert(p.vars.height == x);
		
		double p_plant_survival = exp(-(*istate_ex));

		*irates_ex++ = p.vars.mortality_dt;	// mortality
		*irates_ex++ = p.vars.fecundity_dt * env.patch_survival(t) * p_plant_survival; // viable_seeds
		*irates_ex++ = p.vars.area_heartwood_dt; // heartwood area
		*irates_ex++ = p.vars.mass_heartwood_dt; // heartwood mass

		return irates_ex;
	
	}

	// For output only
	void setState(Solver<PlantModel> *S){
		
		auto iset = S->getIterators_state();

		p.vars.mortality      = *iset.get("mort");
		p.vars.fecundity      = *iset.get("fec");
		p.vars.area_heartwood = *iset.get("heart_area");
		p.vars.mass_heartwood = *iset.get("heart_mass");

	}


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





int main(){
	
	initPlantParameters(plant::par);
	
	//plant::Environment env(1);

	plant::Plant p;
	//p.lma = 0.1978791; // 0.0825;
	//plant::par.r_l   = 198.4545; //39.27 / 0.1978791; // JAI: Should be 39.27/lma;
	//plant::par.k_l = 0.4565855;
	//p.set_height(0.3441948);
	//for (int i=0; i<10000; ++i) p.compute_vars_phys(env);

	cout << p << endl;

	Solver<PlantModel> S(vector<double> (1,p.vars.height), SOLVER_CM);
	S.createSizeStructuredVariables({"mort", "fec", "heart_area", "heart_mass"});

	PlantModel M;
	M.p = p;
	S.setModel(&M);
	//S.createSizeStructuredVariables({"mort", "fec", "heart", "sap"});
	//S.print();
	S.setInputNewbornDensity(M.input_seed_rain);

	S.initialize();
	//M.computeEnv(0, S.state, &S);
	//S.calcRates_CM(1, S.state, S.rates);
	//S.calcRates_extra(1, S.state, S.rates);
	S.print();
	cout << "state: "; for (auto s :S.state) cout << s << " "; cout << "\n";


	//S.addCohort_CM();
	//S.print();
	//cout << "state: "; for (auto s :S.state) cout << s << " "; cout << "\n";



//    double t0 = 0, tf = 50, dt = 1;
//    size_t nsteps = (tf-t0)/dt;
////	vector <double> heights;// = {p.vars.height};

	vector <double> times = generateDefaultCohortSchedule(105.32);
	for (auto t : times) cout << t << " "; cout << endl;

#define OUTPUT_TEXT

#ifdef OUTPUT_TEXT
	ofstream fout("patch_full_hts.txt");
	ofstream fout_ld("patch_full_lds.txt");
	ofstream fout_m("patch_full_m.txt");
	ofstream fout_vs("patch_full_vs.txt");
	ofstream fout_ha("patch_full_ha.txt");
	ofstream fout_hm("patch_full_hm.txt");
	ofstream fli("light_profile_ind_plant.txt");
#endif
	for (size_t i=0; i < times.size(); ++i){

		S.step_to(times[i]);		

#ifdef OUTPUT_TEXT
		cout << times[i] << " " << S.xsize() << " " << M.env.light_profile.npoints << " | " << M.nrc << " " << M.ndc << "\n";
		//S.print();

		vector<double> xl = seq(0, 20, 200);
		for (auto h : xl) fli << M.env.canopy_openness(h) << "\t";
		fli << endl;

		fout << times[i] << "\t";
		fout_ld << times[i] << "\t";
		fout_m << times[i] << "\t";
		fout_vs << times[i] << "\t";
		fout_ha << times[i] << "\t";
		fout_hm << times[i] << "\t";
		//fout_ha << times[i] << "\t";
		auto iset = S.getIterators_state();
		auto& itx = iset.get("X");
		auto& itu = iset.get("u");
		auto& ite = iset.get("mort");
		for (iset.rbegin(); !iset.rend(); --iset){
			fout << *itx << "\t";
			fout_ld << *itu << "\t";
			fout_m << *ite << "\t";
			fout_vs << *next(ite,1) << "\t";
			fout_ha << *next(ite,2) << "\t";
			fout_hm << *next(ite,3) << "\t";
			
		}
		fout << "\n";
		fout_ld << "\n";
		fout_m << "\n";
		fout_vs << "\n";
		fout_ha << "\n";
		fout_hm << "\n";
#endif
	}
	
#ifdef OUTPUT_TEXT
	fli.close();
	fout.close();
	fout_ld.close();
	fout_m.close();
	fout_vs.close();
	fout_ha.close();
	fout_hm.close();
#endif
	cout << "derivative computations requested/done: " << M.nrc << " " << M.ndc << endl;

	auto iset = S.getIterators_state();
	auto& itf = iset.get("fec");
	vector <double> fec_vec;
	fec_vec.reserve(S.xsize());
	iset.rbegin();
	for (int i=0; !iset.rend(); --iset, ++i){
		double patch_age_density = M.env.patch_age_density(times[i]);
		double S_D = 0.25;
		double output_seeds = M.input_seed_rain * S_D * patch_age_density * (*itf);
		cout << times[i] << " " << M.input_seed_rain << " " << S_D << " " << patch_age_density << " " << (*itf) << " | " << output_seeds << endl;
		fec_vec.push_back(output_seeds);
	}
	cout << "Seed rain out = " << pn::integrate_trapezium(times, fec_vec) << endl;
		

}

