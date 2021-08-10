#ifndef DEMO_RED_MODEL_H
#define DEMO_RED_MODEL_H


class LightEnvironment{
	
	double E = 0;

	public:
	double evalEnv(double x, double t){
		return E;
	}

	// This function must do any necessary precomputations to facilitate evalEnv()
	// Therefore, this should calculate env for all X when it is a function of X
	// In such a case, the solver's SubdivisionSpline can be ussed
	// Note: The state vector in the solver will not be updated until the RK step is completed. 
	// Hence, explicitly pass the state to this function.
	template<class Model>
	void computeEnv(double t, vector<double> &state_vec, Solver<Model,LightEnvironment> * S){
		//             _xm 
		// Calculate _/ w(z,t)u(z,t)dz
		//         xb
		auto w = [](double z, double t) -> double {
			return 0.396*pow(z, 0.749)/10000;
		};
		E = S->integrate_x(w, t, state_vec, 0);
	}

};



class REDModel{
	public:

	double input_seed_rain = 1;	


	int nrc = 0; // number of evals of compute_vars_phys() - derivative computations actually done by plant
	int ndc = 0; // number of evals of mortality_rate() - derivative computations requested by solver
	int nbc = 0;

	double a0 = 0.396;
	double phiA = 0.749;
	double g0 = 0.0838;
	double phiG = 0.7134;
	double m0 = 1;
	double mort = 0.035;
	double alpha = 0.1;
	double mu0;

	REDModel() {
		mu0 = mort*m0/g0;
	}

	double initDensity(double x, LightEnvironment * env){
		return 100/pow(x,4);
	}



	double establishmentProbability(double t, LightEnvironment * env){
		return 1;
	}

	double growthRate(double x, double t, LightEnvironment * env){
		++nrc;
		return g0*pow(x,phiG);	
	}

	double mortalityRate(double x, double t, LightEnvironment * env){
		++ndc;
		return mort;
	}

	double birthRate(double x, double t, LightEnvironment * env){
		++nbc;
		return 0.1/0.9*g0*pow(x,phiG)*(1-env->evalEnv(x,t));
	}


	vector<double> initStateExtra(double x, double t, LightEnvironment * env){
		return vector<double>();
	}
		
	vector<double>::iterator calcRates_extra(double x, double t, LightEnvironment * env, 
											vector<double>::iterator istate, vector<double>::iterator irates){
		return vector<double>().begin();
	}

};


#endif
