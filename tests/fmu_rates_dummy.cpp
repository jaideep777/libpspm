#include <vector>
#include <iostream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model.h"

int main(){

	Solver<TestModel> S(25, 0,1, SOLVER_FMU);
	
	TestModel M;
	
	S.setModel(&M);
	S.initialize();
	S.print();	

	M.computeEnv(0,S.state, &S);	
	S.calcRates_FMU(1, S.state, S.rates);  // dummy rates calc rates(X=X0, U=U0, t=1, E=E(U0))
//	S.step_to(1);

	vector <double> rates_exp = {
		-0.355862928, -0.450926505, -0.417535387, -0.354046196, 
		-0.301067280, -0.256580544, -0.219012561, -0.187124217,
	    -0.159930669, -0.136642529, -0.116622139, -0.099350697,
	    -0.084403290, -0.071429752, -0.060139860, -0.050291781,
	    -0.041683001, -0.034143136, -0.027528209, -0.021716055,
	    -0.016602617, -0.012098936, -0.008128703, -0.004910899,
	    -0.001250216};
	for (int i=0; i< S.rates.size(); ++i){
		//cout << S.rates[i] << " " << rates_exp[i] << endl;
		if ( fabs(S.rates[i] - rates_exp[i]) > 1e-4) return 1;
	}
	
	return 0;
	
}


