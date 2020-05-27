#include <vector>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <sys/time.h>
#include <chrono>
#include "../include/plant.h"
#include "../include/cohort.h"
#include "../include/patch.h"
#include "../include/gsl_ode_solver.h"


using namespace std;
using namespace plant;


// -------------------------------------------------------------------------------------------    
//     This code reproduces the growth of metacommunity as with the following code:
//
//     p0 <- scm_base_parameters("FF16")
//     p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FALSE)
//     p1$seed_rain = 200
//     data1 <- run_scm_collect(p1)
//     
//     t <- data1$time
//     h <- data1$species[[1]]["height", , ]
//     t <- data1$time
//     h <- data1$species[[1]]["height", , ]
//     matplot(t, h, lty=1, col=make_transparent("black", 0.25), type="l",
//             las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
// -------------------------------------------------------------------------------------------    



double fixed_canopy_openness(double height){
//  const bool within_canopy = height <= light_environment.max();
//  return within_canopy ? light_environment.eval(height) : 1.0;
	return 1;
}


int main(){

	initPlantParameters(par);
	
	Environment env(1);
	env.light_profile.print();
	cout << "~~~~~~~" << endl;
	
	Patch P(200, env);
	P.compute_light_profile(env);
	cout << P << endl;

	int rate_eval_count = 0;
	auto func_lambda = [&P, &env, &rate_eval_count](double t, const double y[], double f[]) -> int {
		env.time = t;
		P.set_state(y);
		P.compute_vars_phys(env);
		P.get_rates(f);
								   
		++rate_eval_count;
		return GSL_SUCCESS;
	};


	ODE_Solver<decltype(func_lambda)> sys(func_lambda, 6);

	vector <double> state0 = P.get_state();

	sys.initialize(state0);

	double t0 = 0, tf = 50, dt = 1;
	size_t nsteps = (tf-t0)/dt;


	ofstream fout("patch_full_hts.txt");
	ofstream fout_ld("patch_full_lds.txt");
	ofstream fout_vs("patch_full_vs.txt");
//	ofstream fout_le_z("le_z.txt");		
//	ofstream fout_le_co("le_co.txt");		
	ofstream fout_m("patch_full_ms.txt");
	ofstream fout_ha("patch_full_ha.txt");
	ofstream fout_hm("patch_full_hm.txt");

	cout << "Step " << 0 << ", t = " << P.cohort_introduction_times[0] << ", " << P.cohorts.size() << " cohorts." <<  endl;
	fout << env.time << "\t";
	fout_ld << env.time << "\t";
	fout_vs << env.time << "\t";
	fout_m << env.time << "\t";
	fout_ha << env.time << "\t";
	fout_hm << env.time << "\t";
	for (auto &c : P.cohorts){
		Plant &p = c.plant;
		
		fout << 
			p.vars.height         << "\t";// <<

		fout_ld << 
			c.log_density         << "\t";// <<

		fout_vs << 
			c.viable_seeds         << "\t";// <<

		fout_m << 
			p.vars.mortality         << "\t";// <<

		fout_ha << 
			p.vars.area_heartwood         << "\t";// <<

		fout_hm << 
			p.vars.mass_heartwood         << "\t";// <<

	}						
	fout << endl;
	fout_ld << endl;
	fout_vs << endl;
	fout_m << endl;
	fout_ha << endl;
	fout_hm << endl;




	for (size_t i=1; i < P.cohort_introduction_times.size(); ++i){


//		cout << "Light Env Profile Calc: ";
//		auto t1 = std::chrono::steady_clock::now();
//		for (int i=0; i<50; ++i){
//			double z = i*P.height_max()/49.0;
//			double co = env.canopy_openness(z);
//			fout_le_z << setprecision(8) << z << "\t";
//			fout_le_co << setprecision(8) << co << "\t";
//		}
//		fout_le_z << endl;
//		fout_le_co << endl;
//		auto t2 = std::chrono::steady_clock::now();
//		cout << "[" << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/1e3 << " usec]" << std::endl;


		sys.step_to(P.cohort_introduction_times[i]);
		//cout << "Required rate evaluations: " << rate_eval_count << endl;
		rate_eval_count = 0;

//		cout << "Patch seed rain = " << P.seed_rain(env) << endl;
		
		env.time = P.cohort_introduction_times[i];

//		if (i<40){
		P.introduce_cohort(env);

		vector <double> new_state = P.get_state();
		sys.resize(P.cohorts.size()*6);
		sys.initialize(new_state, env.time);
//		}

//		cout << "Step " << i << ", t = " << P.cohort_introduction_times[i] << ", " << P.cohorts.size() << " cohorts." <<  endl;
		fout << env.time << "\t";
		fout_ld << env.time << "\t";
		fout_vs << env.time << "\t";
		fout_m << env.time << "\t";
		fout_ha << env.time << "\t";
		fout_hm << env.time << "\t";
		for (auto &c : P.cohorts){
			Plant &p = c.plant;
			
			fout << 
				p.vars.height         << "\t";// <<

			fout_ld << 
				c.log_density         << "\t";// <<

			fout_vs << 
				c.viable_seeds         << "\t";// <<

			fout_m << 
				p.vars.mortality         << "\t";// <<

			fout_ha << 
				p.vars.area_heartwood         << "\t";// <<

			fout_hm << 
				p.vars.mass_heartwood         << "\t";// <<

		}						
		fout << endl;
		fout_ld << endl;
		fout_vs << endl;
		fout_m << endl;
		fout_ha << endl;
		fout_hm << endl;

	
	}


//	env.light_profile.printToFile("light_profile.txt");

	fout.close();
	fout_ld.close();
	fout_vs.close();


	cout << "Seed rain: " << P.input_seed_rain << " ---> " << P.seed_rain(env) << endl;


	return 0;
}


