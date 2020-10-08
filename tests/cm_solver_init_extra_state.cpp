#include <vector>
#include <iostream>
#include <cmath>
#include "solver.h"
#include <fstream>
using namespace std;

#include "test_model_2.h"


int main(){

	Solver<TestModel> S(10, 0,1, SOLVER_CM);
	
	TestModel M;
	
	S.setModel(&M);
	S.createSizeStructuredVariables({"mort", "vs", "heart", "sap"});

	S.initialize();
	M.computeEnv(0, S.state, &S);
	S.calcRates_CM(1, S.state, S.rates);
	S.calcRates_extra(1, S.state, S.rates);
	S.print();

	ifstream fin("tests/test_data/cm_extra_state_rates_J10_testmodel2.txt");
	for (int i=0; i<S.state.size(); ++i){
		double d;
		fin >> d;
		//cout << d << " " << S.state[i] << endl;
		if (abs(d - S.state[i]) > 1e-4 && !isinf(S.state[i])){
			return 1;
		}
	}
	for (int i=0; i<S.rates.size(); ++i){
		double d;
		fin >> d;
		if (abs(d - S.rates[i]) > 1e-4){
			//cout << d << " " << S.rates[i] << endl;
			return 1;
		}
	}
	fin.close();

	vector<double> xbreaks = {0, 0.25, 0.5, 0.75, 1};
	S.resetState(xbreaks);
	S.initialize();
	S.print();


	
	return 0;
}


