#include <vector>
#include <iostream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model.h"

int main(){

	Solver<TestModel> S(25, 0,1, SOLVER_CM);
	
	TestModel M;
	
	S.setModel(&M);
	S.initialize();
	//S.print();	
	
	M.computeEnv(0,&S);

	if (fabs(M.evalEnv(0,0) - 0.3821924) > 1e-6) return 1;
	
	return 0;
	
}

