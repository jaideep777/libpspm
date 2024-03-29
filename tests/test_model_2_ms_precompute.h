#ifndef __PSPM_TEST_TEST_TEST_MODEL_H_
#define __PSPM_TEST_TEST_TEST_MODEL_H_



class Environment : public EnvironmentBase{
	
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
	void computeEnv(double t, Solver * S, std::vector<double>::iterator s, std::vector<double>::iterator dsdt){
		//             _xm 
		// Calculate _/ w(z,t)u(z,t)dz
		//         xb
		auto w = [S](int i, double t) -> double {
			double z = S->species_vec[0]->getX(i);
			if (z <= 1.0/3) 
				return 1;
			else if (z > 1.0/3 && z <= 2.0/3) 
				return pow(2-3*z, 3)*(54*z*z-27*z+4);
			else 
				return 0;
		};
		E = S->integrate_x(w, t, 0);
	}

};



class Plant{
	public:
	double height;
	double mortality;
	double viable_seeds;
	double heart_mass;
	double sap_mass;

	std::vector<std::string> varnames = {"mort", "vs", "heart", "sap"};

	Plant(double h){
		height = h;
	}

	std::vector<double> calcRates(){
		std::vector<double> rates(4);
		rates[0] = -2;
		rates[1] = 0;
		rates[2] = -20*height;
		rates[3] = -30*height;
		return rates;
	}
};





class TestModel : public Plant{
	public:
	double sc = 10;

	// precomputed demographic rates
	double g,m,f,se;
	
	TestModel() : Plant(0) {}

	double init_density(double x, void * _env, double bf){
		return pow(1-x,2)/pow(1+x,4) + (1-x)/pow(1+x,3);
	}
	
	void preCompute(double x, double t, void * _env){
		Environment* env = (Environment*)_env;
	
		// growth rate
		double E = env->evalEnv(x,t);
		double a = 0.16+0.22*exp(-0.225*t*t);
		g = 0.225*(1-x*x)*(E/(1+E*E))*t*(1+a*a)/a;

		m = 1.35*t*E/a;
		
		double oneplusa = 1+a;
		double n1 = 0.225*t*x*x*(1-x)*(1-x)*E/(1+E)/(1+E)*oneplusa*oneplusa/a;
		double n2 = (1+exp(-0.225*t*t))/(61-88*log(2)+(38*log(2)-79.0/3)*exp(-0.225*t*t));
		f = n1*n2;

		se = 1;

	}

	double growthRate(double x, double t, void * _env){
		return g;
	}

	double mortalityRate(double x, double t, void * _env){
		return m;
	}

	double birthRate(double x, double t, void * _env){
		return f;
	}

	double establishmentProbability(double t, void  * _env){
		return se;
	}


	void set_size(double _x){
		height = _x;
		mortality = 0.1*height + 1e-12; //exp(-height);	
		viable_seeds = 100*height + 1e-13;
		heart_mass = 1000*height + 1e-14;
		sap_mass = 10*height + 1e-15;
	}

	void init_state(double t, void * env){
	}

	std::vector<double>::iterator set_state(std::vector<double>::iterator &it){
		mortality    = *it++;
		viable_seeds = *it++;
		heart_mass   = *it++;
		sap_mass     = *it++;
		return it;
	}
	std::vector<double>::iterator get_state(std::vector<double>::iterator &it){
		*it++ = mortality;
		*it++ = viable_seeds;
		*it++ = heart_mass;
		*it++ = sap_mass;
		return it;
	}

	std::vector<double>::iterator get_rates(std::vector<double>::iterator &it){
		std::vector<double> r = calcRates();
		*it++ = r[0];
		*it++ = r[1];
		*it++ = r[2];
		*it++ = r[3];
		return it;
	}

	void print(ostream &out = cout){
		out << mortality << "\t" << viable_seeds << "\t" << heart_mass << "\t" << sap_mass << "\t";
	}


};


#endif
