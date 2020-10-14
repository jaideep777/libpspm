#ifndef __PSPM_TEST_TEST_TEST_MODEL_H_
#define __PSPM_TEST_TEST_TEST_MODEL_H_

class Environment{
	
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
	void computeEnv(double t, vector<double> &state_vec, Solver<Model,Environment> * S){
		//            _xm 
		// Calculate / w(z,t)u(z,t)dz
		//        xb`
		auto w = [](double z, double t) -> double {
			if (z <= 1.0/3) 
				return 1;
			else if (z > 1.0/3 && z <= 2.0/3) 
				return pow(2-3*z, 3)*(54*z*z-27*z+4);
			else 
				return 0;
		};
		E = S->integrate_x(w, t, state_vec, 0);
	}

};



class Plant{
	public:
	double height;
	double mortality;
	double viable_seeds;
	double heart_mass;
	double sap_mass;
	
	Plant(double h){
		height = h;
	}

	void init_state(){
		mortality = 3*exp(-height);	
		viable_seeds = 100*height;
		heart_mass = 1000*height*height*height;
		sap_mass = 10*sqrt(height);
	}

	vector<double> calcRates(){
		vector<double> rates(4);
		rates[0] = -2;
		rates[1] = 0;
		rates[2] = -20*height;
		rates[3] = -30*height;
		return rates;
	}
};





class TestModel{
	public:
	double sc = 10;
	
	double initDensity(double x){
		return pow(1-x,2)/pow(1+x,4) + (1-x)/pow(1+x,3);
	}
	
	double growthRate(double x, double t, Environment* env){
		double E = env->evalEnv(x,t);
		double a = 0.16+0.22*exp(-0.225*t*t);
		return 0.225*(1-x*x)*(E/(1+E*E))*t*(1+a*a)/a;
	}

	double mortalityRate(double x, double t, Environment* env){
		double E = env->evalEnv(x,t);
		double a = 0.16+0.22*exp(-0.225*t*t);
		return 1.35*t*E/a;
	}

	double birthRate(double x, double t, Environment* env){
		double E = env->evalEnv(x,t);
		double oneplusa = 1.16+0.22*exp(-0.225*t*t);
		double a = 0.16+0.22*exp(-0.225*t*t);
		double n1 = 0.225*t*x*x*(1-x)*(1-x)*E/(1+E)/(1+E)*oneplusa*oneplusa/a;
		double n2 = (1+exp(-0.225*t*t))/(61-88*log(2)+(38*log(2)-79.0/3)*exp(-0.225*t*t));
		return n1*n2;
	}

	double establishmentProbability(double t){
		return 1;
	}
	
	// optional functions, if extra size-structured variables are desired
	vector<double> initStateExtra(double x, double t){
		Plant p(x);
		p.init_state();
		vector<double> sv;
	    sv.reserve(4);	
		sv.push_back(p.mortality); 
		sv.push_back(p.viable_seeds); 
		sv.push_back(p.heart_mass); 
		sv.push_back(p.sap_mass);
		return sv;
	}

	vector<double>::iterator calcRates_extra(double t, double x, 
											vector<double>::iterator istate, vector<double>::iterator irates){
		
		Plant p(x);
		p.init_state();
		vector<double> r = p.calcRates();
		*irates++ = r[0];
		*irates++ = r[1];
		*irates++ = r[2];
		*irates++ = r[3];

		return irates;
	
	}

	//void computeEnv(double t, vector<double> &state_vec, Solver<TestModel> * S){
	//}
};


#endif
