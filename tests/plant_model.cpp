#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
#include "solver.h"
#include "plant/plant.h"


class PlantModel{
	public:
	
	plant::Plant p;
	plant::Environment env;

	int nrc = 0; // number of evals of compute_vars_phys() - derivative computations actually done by plant
	int ndc = 0; // number of evals of mortality_rate() - derivative computations requested by solver

	PlantModel() : env(1) {
	
	}

	double initDensity(double x){
		return 1;
	}

	double evalEnv(double x, double t){
		return 1;
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
	void computeEnv(double t, vector<double> &state_vec, Solver<PlantModel> * S){
		//            _xm 
		// Calculate / w(z,t)u(z,t)dz
		//        xb`
	}

	double growthRate(double x, double t){
		if (p.vars.height != x){
			p.set_height(x);
			p.compute_vars_phys(env);
			++nrc;
		}
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
	vector<double> initStateExtra(double x){
		vector<double> sv;
		sv.reserve(4);	
		sv.push_back(0); 
		sv.push_back(0); 
		sv.push_back(0); 
		sv.push_back(0);
		return sv;
	}


	vector<double>::iterator calcRates_extra(double t, double x, vector<double>::iterator irates_ex){
		
		assert(p.vars.height == x);
		
		*irates_ex++ = p.vars.mortality_dt;
		*irates_ex++ = p.vars.fecundity_dt;
		*irates_ex++ = p.vars.area_heartwood_dt;
		*irates_ex++ = p.vars.mass_heartwood_dt;

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

int main(){
	
	initPlantParameters(plant::par);
	
	plant::Environment env(1);

	plant::Plant p;
	p.lma = 0.1978791;
	plant::par.r_l   = 198.4545; //39.27 / 0.1978791; // JAI: Should be 39.27/lma;
	plant::par.k_l = 0.4565855;

	p.set_height(0.3441948);
	for (int i=0; i<10000; ++i) p.compute_vars_phys(env);

	cout << p << endl;

	Solver<PlantModel> S(vector<double> (1,p.vars.height), SOLVER_CM);
	S.createSizeStructuredVariables({"mort", "fec", "heart_area", "heart_mass"});

	PlantModel M;
	S.setModel(&M);
	//S.createSizeStructuredVariables({"mort", "fec", "heart", "sap"});
	S.print();

	S.initialize();
	M.computeEnv(0, S.state, &S);
	S.calcRates_CM(1, S.state, S.rates);
	S.calcRates_extra(1, S.state, S.rates);
	S.print();




	double t0 = 0, tf = 50, dt = 1;
	size_t nsteps = (tf-t0)/dt;
//	vector <double> heights;// = {p.vars.height};

	ofstream fout("ind_plant.txt");

	for (size_t i=0; i < nsteps; ++i){

		S.step_to(i*dt);		
		
		M.setState(&S);

		fout << i*dt << "\t" <<
			M.p.vars.height         << "\t" <<
			M.p.vars.mortality      << "\t" <<
			M.p.vars.fecundity      << "\t" <<
			M.p.vars.area_heartwood << "\t" <<
			M.p.vars.mass_heartwood << endl;
						
//		cout << p.vars.fecundity << " ";		
//		heights.push_back(p.vars.height);
	}
	
	cout << M.p << endl;
	cout << "derivative computations requested/done: " << M.nrc << " " << M.ndc << endl;


}

