#include <vector>
#include <iostream>
#include <cmath>
#include "solver.h"
using namespace std;


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

};



class Model{
	public:
	double sc = 10;
	
	double calcIC(double x){
		return sc*x;
	}

	vector<double> initStateExtra(double x){
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

};


int main(){

	Solver<Model> S(10, 0,1, SOLVER_CM);
	
	Model M;
	
	S.setModel(&M);
	S.createSizeStructuredVariables({"mort", "vs", "heart", "sap"});

	S.initialize();
	S.print();	

	const double * X = S.getX();
	for (int i=0; i<S.xsize(); ++i){
		cout << X[i] << " ";
	}
	cout << endl;
	
	vector<double> xbreaks = {0.4, 0.8};
	S.resetState(xbreaks);
	S.initialize();
	S.print();


	
	return 0;
}


