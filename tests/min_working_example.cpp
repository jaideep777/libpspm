#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include <individual_base.h>
#include <solver.h>

using namespace std;

class Environment : public EnvironmentBase{
	public:
	double E = 0;

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



class Individual : public IndividualBase{
	public:

	double init_density(double x, void * _env, double bf){
		return pow(1-x,2)/pow(1+x,4) + (1-x)/pow(1+x,3);
	}

	double growthRate(double x, double t, void * _env){
		Environment* env = (Environment*)_env;
		double E = env->E;
		double a = 0.16+0.22*exp(-0.225*t*t);
		return 0.225*(1-x*x)*(E/(1+E*E))*t*(1+a*a)/a;
	}

	double mortalityRate(double x, double t, void * _env){
		Environment* env = (Environment*)_env;
		double E = env->E;
		double a = 0.16+0.22*exp(-0.225*t*t);
		return 1.35*t*E/a;
	}

	double birthRate(double x, double t, void * _env){
		Environment* env = (Environment*)_env;
		double E = env->E;
		double oneplusa = 1.16+0.22*exp(-0.225*t*t);
		double a = 0.16+0.22*exp(-0.225*t*t);
		double n1 = 0.225*t*x*x*(1-x)*(1-x)*E/(1+E)/(1+E)*oneplusa*oneplusa/a;
		double n2 = (1+exp(-0.225*t*t))/(61-88*log(2)+(38*log(2)-79.0/3)*exp(-0.225*t*t));
		return n1*n2;
	}

};


int main(){
	Environment E;
	Species<Individual> spp;

	Solver S(SOLVER_FMU, "rk45ck");

	S.addSpecies(25, 0, 1, false, &spp, 0, -1);
	S.setEnvironment(&E);

	S.resetState();
	S.initialize();

	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		cout << S.current_time << " " << S.u0_out(t)[0] << "\n";
	}

}


