#include <vector>
#include <iostream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model_2_ms.h"

int main(){

	TestModel M;
	
	Solver<TestModel,Environment> S(SOLVER_CM);
	S.use_log_densities = false;
	S.control.cm_grad_dx = 0.001;
	S.addSpecies(25, 0, 1, false, &M);
	S.resetState();
	S.initialize();
	//S.print();
	
	Environment E;
	E.computeEnv(0,S.state,&S);
	//cout << E.evalEnv(0,0) << endl;

	S.setEnvironment(&E);
	S.calcRates_CM(1, S.state, S.rates);  // dummy rates calc rates(X=X0, U=U0, t=1, E=E(U0))
	//S.print();
	//	S.step_to(1);

	vector <double> rates_exp = {
		  0.24871666,  0.24831871,  0.24712487,  0.24513514,  0.24234951,  0.23876799,  0.23439058,  0.22921727,
		  0.22324807,  0.21648298,  0.20892199,  0.20056511,  0.19141234,  0.18146367,  0.17071911,  0.15917866,
		  0.14684231,  0.13371007,  0.11978194,  0.10505792,  0.08953800,  0.07322218,  0.05611048,  0.03820288,
		  0.01949939,  0.00000000, -3.07367787, -2.48964104, -2.02468490, -1.65220631, -1.35210765, -1.10906811,
		 -0.91130994, -0.74970832, -0.61714231, -0.50801705, -0.41790887, -0.34329937, -0.28137461, -0.22987235,
		 -0.18696499, -0.15116931, -0.12127634, -0.09629670, -0.07541760, -0.05796891, -0.04339629, -0.03123970,
		 -0.02111622, -0.01270624, -0.00574230, 0.00000000 };
	

	for (int i=0; i< S.rates.size(); ++i){
		//cout << S.rates[i] << " " << rates_exp[i] << endl;
		if ( fabs(S.rates[i] - rates_exp[i]) > 1e-6) return 1;
	}
	
	return 0;
	
}


