#include <iostream>
#include <solver.h>

#include "test_model_2_ms.h"

int main(){
	
	TestModel M;
	Environment E;

	Solver<TestModel, Environment> S(SOLVER_CM);
	S.setEnvironment(&E);
	S.addSpecies(10, 0, 1, false, &M, {"mort", "fec"}, 10);
	S.addSpecies(5, 0, 0.5, false, &M, {"ha"}, 5);
	S.addSpecies(2, 0, 1, false, &M, {}, 2);
	S.resetState();
	S.initialize();
	E.computeEnv(0, S.state, &S);
	S.print();

	S.addCohort_CM();
	S.print();

	S.addCohort_CM();
	S.print();

	S.removeCohort_CM();
	S.print();

	return 0;
}

