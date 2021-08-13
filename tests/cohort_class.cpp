#include <iostream>
#include <cohort.h>
#include <species.h>
#include <solver.h>
#include <vector>
using namespace std;


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


class Environment{
	double E = 0;
};

int main(){
	Cohort<Plant> C;
	
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

	std::cout << " " << C.lma << " " << C.x << " " << C.u << " " << C.id << "\n";

	Solver<Environment> sol(SOLVER_FMU);
	sol.addSpecies(vector<double> {1,2,3,4,5}, &S, 1, 1);
	sol.addSpecies(vector<double> {1.1,2.1,3.1}, &I, 0, 1);
	sol.resetState();
	sol.initialize();
	sol.print();	

}



