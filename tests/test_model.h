#ifndef __PSPM_TEST_TEST_TEST_MODEL_H_
#define __PSPM_TEST_TEST_TEST_MODEL_H_

class TestModel{
	public:
	double sc = 10;
	double env;
	
	double calcIC(double x){
		return pow(1-x,2)/pow(1+x,4) + (1-x)/pow(1+x,3);
	}

	vector<double> initStateExtra(double x){
	
	}

	double evalEnv(double x, double t){
		return env;
	}
	
	// This function must do any necessary precomputations to facilitate evalEnv()
	// Therefore, this should calculate env for all X when it is a function of X
	// In such a case, the solver's SubdivisionSpline can be ussed
	// Note: The state vector in the solver will not be updated until the RK step is completed. 
	// Hence, explicitly pass the state to this function.
	void computeEnv(double t, vector<double> &state_vec, Solver<TestModel> * S){
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
		env = S->integrate_x(w, t, state_vec, 1);
	}

	double growthRate(double x, double t){
		double E = evalEnv(x,t);
		double a = 0.16+0.22*exp(-0.225*t*t);
		return 0.225*(1-x*x)*(E/(1+E*E))*t*(1+a*a)/a;
	}

	double mortalityRate(double x, double t){
		double E = evalEnv(x,t);
		double a = 0.16+0.22*exp(-0.225*t*t);
		return 1.35*t*E/a;
	}

	double birthRate(double x, double t){
		double E = evalEnv(x,t);
		double oneplusa = 1.16+0.22*exp(-0.225*t*t);
		double a = 0.16+0.22*exp(-0.225*t*t);
		double n1 = 0.225*t*x*x*(1-x)*(1-x)*E/(1+E)/(1+E)*oneplusa*oneplusa/a;
		double n2 = (1+exp(-0.225*t*t))/(61-88*log(2)+(38*log(2)-79.0/3)*exp(-0.225*t*t));
		return n1*n2;
	}
		
};



#endif
