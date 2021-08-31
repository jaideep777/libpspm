#include <iostream>
#include <cohort.h>
#include <species.h>
#include <solver.h>
#include <vector>
#include <environment_base.h>
using namespace std;


class LightEnv : public EnvironmentBase{
	public:
	double E = 0.95;

	void computeEnv(double t, Solver * S){
		E = 0.95 + 0.05*t;
	}

};



class Plant {
	public:
	double lma = 30;

	double height = -99;
	double crown_area = -99;
	double root_mass = -99;

	vector<string> varnames = {"lma|", "ht", "cr", "root"};

	Plant(){
		lma = 10;
	}

	Plant(double _lma){
		lma = _lma;
	}
	
	void set_size(double x){
		height = x;
		crown_area = x*x;
	}

	void preCompute(double x, double t, void * env){
	}

	double growthRate(double x, double t, void * env){
		cout << "in g: " << x << " " << t << " " << ((LightEnv*)env)->E << "\n";
		return x*((LightEnv*)env)->E;
	}
	double mortalityRate(double x, double t, void * env){
		return -0.5;
	}
	double birthRate(double x, double t, void * env){
		return 1;
	}
	
	double establishmentProbability(double t, void  * _env){
		return 1;
	}

	double init_density(double x, void * _env){
		return 5/(x+0.5);
	}

	void init_state(double t, void * env){
		root_mass = 10;
	}

	vector<double>::iterator set_state(vector<double>::iterator &it){
		root_mass = *it++;
		return it;
	}
	vector<double>::iterator get_state(vector<double>::iterator &it){
		*it++ = root_mass;
		return it;
	}

	vector<double>::iterator get_rates(vector<double>::iterator &it){
		*it++ = -0.1*(root_mass-1);
		return it;
	}

	void print(std::ostream& out = std::cout){
		out << lma << "|\t" << height << "\t" << crown_area << "\t" << root_mass << "\t";
	}
};


class Insect {
	public:
	double wingspan = 10;
	
	double g = -1, m = -1, f = -1;

	vector<string> varnames = {"wing", "f"};
	
	double init_density(double x, void * _env){
		return x/10;
	}

	void init_state(double t, void * env){
	}

	void set_size(double x){
		wingspan = 2*x;	
	}

	void preCompute(double x, double t, void * env){
		f = x*100+t;
		g = x*50*t;
		m = 0.1;
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
	
	double growthRate(double x, double t, void * env){
		return g;
	}
	double mortalityRate(double x, double t, void * env){
		return m;
	}
	double birthRate(double x, double t, void * env){
		return f;
	}

	double establishmentProbability(double t, void  * _env){
		return 1;
	}

	void print(std::ostream& out = std::cout){
		out << wingspan << "\t" << f << "\t";
	}
};



int main(){

	Plant P;
	cout << "P: "; P.print(); cout << "\n";

	P.lma = 70;
	Cohort<Plant> C1;
	C1 = Cohort<Plant> (P);
	cout << "C1: "; C1.print(); cout << "\n";


	Species<Plant> Sv(vector<double> {1,2,3,4,5});
	Sv.setX(0, 1.5);
	Sv.setX(2, 3.5);
	Sv.print();

	Species<Plant> S(P);  // create species S with prototype P

	Species<Insect> I(vector<double> {1.1,2.1,3.1});

	Species_Base * S1 = &S;
	//S1->print();
	cout << "S size: " << S1->get_maxSize() << "\n";

	Species_Base * S2 = &I;
	//S2->print();
	cout << "I size: " << S2->get_maxSize() << "\n";


	Solver sol(SOLVER_EBT);
	sol.addSpecies(vector<double> {1,2,3,4,5}, &S, 1, 1);
	sol.addSpecies(vector<double> {1.1,2.1,3.1}, &I, 0, 1);
	sol.resetState();
	sol.initialize();
	//sol.addCohort_EBT();
	sol.print();	

	I.setX(2, 0.05);
	I.setU(2, 0.2);
	sol.print();	
	cout << "Insect BF = " << sol.calcSpeciesBirthFlux(1,0) << "\n";
	
	sol.preComputeSpecies(1,0);
	sol.print();	
	cout << "Insect BF (after precompute) = " << sol.calcSpeciesBirthFlux(1,0) << "\n";

	LightEnv env;
	sol.setEnvironment(&env);

	Cohort<Plant> C;
	C.print(); cout << "\n";
	C.set_size(10);
	C.print(); cout << "\n";
	cout << C.growthRate(C.height, 0, sol.env) << "\n";

}



