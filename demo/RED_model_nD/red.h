#ifndef DEMO_RED_MODEL_H
#define DEMO_RED_MODEL_H

#include <individual_base.h>

class LightEnvironment : public EnvironmentBase{
	
	double E = 0;

	public:
	double evalEnv(double xn, double t){
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
			std::vector<double> z = S->species_vec[0]->getXn(i);
			return 0.396*pow(z[0], 0.749)/10000;
		};
		E = S->integrate_x(w, t, 0);
	}

};



class RED_Plant : public IndividualBase{
	public:

	//double input_seed_rain = 1;	
	vector <string> varnames;

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
	double beta = 0.002;
	double mu0;

	RED_Plant() {
		mu0 = mort*m0/g0;
	}

	void set_size(double _x){
	}

	void set_size(std::vector<double> _x){
	}

	double init_density(double x, void * env, double input_seed_rain){
		return 100/pow(x,4);
	}

	double init_density(std::vector<double> xn, void * env, double input_seed_rain){
		return 100/pow(xn[0],4);
	}

	void preCompute(double x, double t, void * env){
	}

	void preCompute(std::vector<double> xn, double t, void * env){
	}

	double establishmentProbability(double t, void * env){
		return 1;
	}

	double growthRate(double x, double t, void * env){
		++nrc;
		return g0*pow(x,phiG);	
	}

	std::vector<double> growthRate(std::vector<double> xn, double t, void * env){
		++nrc;
		// std::cout << "growth Rate Gradient compute" << std::endl;
		// std::cout << xn << std::endl;
		std::vector<double> vec_out;
		double x0 = g0*pow(xn[0],phiG);
		// std::cout << x0 << std::endl;	
		vec_out.push_back(x0);
		double x1 = -beta * xn[1];
		// std::cout << x1 << std::endl;
		vec_out.push_back(x1);
		// std::cout <<x0 << "   " <<  x1 << std::endl;
		// std::cout << "growth Rate Gradient compute: END" << std::endl;

		return {x0, x1};
	}

	double mortalityRate(double x, double t, void * env){
		++ndc;
		return mort;
	}

	double mortalityRate(std::vector<double> xn, double t, void * env){
		++ndc;
		return mort;
	}

	double birthRate(double x, double t, void * env){
		++nbc;
		if (x < 0) throw std::runtime_error("x is negative");
		LightEnvironment* env1 = (LightEnvironment*)env;
		return 0.1/0.9*g0*pow(x,phiG)*(1-env1->evalEnv(x,t));
	}

	double birthRate(std::vector<double> xn, double t, void * env){
		++nbc;
		if (xn[0] < 0) throw std::runtime_error("x_0 is negative");
		LightEnvironment* env1 = (LightEnvironment*)env;
		return 0.1/0.9*g0*pow(xn[0],phiG)*(1-env1->evalEnv(xn[0],t));
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
