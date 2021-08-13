#include <iostream>
#include <cohort.h>
#include <species.h>
#include <solver.h>
#include <vector>
#include <environment.h>
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

	double height;
	double crown_area;
	double root_mass;

	vector<string> varnames = {"ht", "cr", "root"};

	double set_height(double x){
		height = x;
		crown_area = x*x;
		root_mass = 10*x;
	}

	
	double growthRate(double x, double t, void * env){
		cout << "in g: " << x << " " << t << " " << ((LightEnv*)env)->E << "\n";
		return x*((LightEnv*)env)->E;
	}


	double init_density(double x){
		return 5/(x+0.5);
	}

	vector<double>::iterator init_state(double x, vector<double>::iterator &it){
		root_mass = 10*x;
		*it++ = root_mass;
		return it;
	}

	vector<double>::iterator set_state(vector<double>::iterator &it){
		root_mass = *it++;
		return it;
	}

	vector<double>::iterator get_rates(vector<double>::iterator &it){
		*it++ = -0.1*(root_mass-1);
		return it;
	}

	void print(){
		cout << height << "\t" << crown_area << "\t" << root_mass << "\t";
	}
};


class Insect {
	public:
	double wingspan = 10;

	vector<string> varnames = {};
	
	double init_density(double x){
		return x/10;
	}

	vector<double>::iterator init_state(double x, vector<double>::iterator &it){
		return it;
	}

	vector<double>::iterator set_state(vector<double>::iterator &it){
		return it;
	}

	vector<double>::iterator get_rates(vector<double>::iterator &it){
		return it;
	}

	void print(){
	
	}
};



int main(){

	
	Species<Plant> S(vector<double> {1,2,3,4,5});
	S.setX(0, 1.5);
	S.setX(2, 3.5);
	S.print();

	Species<Insect> I(vector<double> {1.1,2.1,3.1});

	Species_Base * S1 = &S;
	S1->print();
	cout << "S size: " << S1->get_maxSize() << "\n";

	Species_Base * S2 = &I;
	S2->print();
	cout << "I size: " << S2->get_maxSize() << "\n";


	Solver sol(SOLVER_FMU);
	sol.addSpecies(vector<double> {1,2,3,4,5}, &S, 1, 1);
	sol.addSpecies(vector<double> {1.1,2.1,3.1}, &I, 0, 1);
	sol.resetState();
	sol.initialize();
	sol.print();	
	
	LightEnv env;
	sol.setEnvironment(&env);

	Cohort<Plant> C;
	C.set_height(10);
	C.print(); cout << "\n";
	cout << C.growthRate(C.height, 0, sol.env) << "\n";

}



