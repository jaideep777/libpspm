#ifndef DEMO_DAPHNIA_MODEL_H
#define DEMO_DAPHNIA_MODEL_H


class Environment : public EnvironmentBase{
	
	public:
	double r = 0.5;
	double K = 3;
	double S = 0;

	double dSdt;
	
	public:
//	double evalEnv(double x, double t){
//		return E;
//	}

//	void calcRatesSystem(double t, vector<double>::iterator S, vector<double>::iterator dSdt){

//	}
	// This function must do any necessary precomputations to facilitate evalEnv()
	// Therefore, this should calculate env for all X when it is a function of X
	// In such a case, the solver's SubdivisionSpline can be ussed
	// Note: The state vector in the solver will not be updated until the RK step is completed. 
	// Hence, explicitly pass the state to this function.
	void computeEnv(double t, Solver * sol, vector<double>::iterator _S, vector<double>::iterator _dSdt){
		//             _xm 
		// Calculate _/ w(z,t)u(z,t)dz
		//         xb
		S = *_S;
		
		auto w = [sol](int i, double t) -> double {
			double z = sol->species_vec[0]->getX(i);
			return z*z;
		};
		double E = sol->integrate_x(w, t, 0);
//		cout << "E = " << E << "\n";
		
		*_dSdt = r*S*(1-S/K) - S/(1+S)*E;
	}

};



class Daphnia{
	public:

	//double input_seed_rain = 1;	
	vector <string> varnames;

	int nrc = 0; // number of evals of compute_vars_phys() - derivative computations actually done by plant
	int ndc = 0; // number of evals of mortality_rate() - derivative computations requested by solver
	int nbc = 0;

	double a = 0.75;
	double mu0 = 0.1;
	
	double xb = 0, xm = 1; 

	Daphnia() {
	}


	void set_size(double _x){
	}

	double init_density(double x, void * env, double input_seed_rain){
		return exp(-8*pow((x-xb)/(xm-xb),3));
	}

	void preCompute(double x, double t, void * env){
	}

	double establishmentProbability(double t, void * env){
		return 1;
	}

	double growthRate(double x, double t, void * env){
		++nrc;
		double S = ((Environment*)env)->S;
		return fmax(S/(1+S) - x, 0);	
	}

	double mortalityRate(double x, double t, void * env){
		++ndc;
		return mu0;
	}

	double birthRate(double x, double t, void * env){
		++nbc;
		double S = ((Environment*)env)->S;
		return a*x*x*S/(1+S);
	}

	void init_state(double t, void * env){
	}
	vector<double>::iterator set_state(vector<double>::iterator &it){
		return it;
	}
	vector<double>::iterator get_state(vector<double>::iterator &it){
		return it;
	}
	vector<double>::iterator get_rates(vector<double>::iterator &it){
		return it;
	}

	void print(std::ostream& out = std::cout){
	}
};


#endif
