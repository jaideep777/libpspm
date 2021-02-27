#include <vector>
#include <iostream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model_2_ms.h"

int main(){

	TestModel M;
	
	Solver<TestModel,Environment> S(SOLVER_CM);
	S.addSpecies(25, 0, 1, false, &M);
	S.resetState();
	S.initialize();
	S.print();
	
	Environment E;
	E.computeEnv(0,S.state,&S);
	cout << E.evalEnv(0,0) << endl;

	if (fabs(E.evalEnv(0,0) - 0.3821924) > 1e-6) return 1;
	
	return 0;
	
}

