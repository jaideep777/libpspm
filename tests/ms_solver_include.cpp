#include <iostream>
#include <solver.h>

#include "test_model_2.h"

int main(){
	
	TestModel M;

	Solver<TestModel> S(SOLVER_CM);
	S.addSpecies(10, 0, 1, false, &M, {"mort", "fec"}, 10);
	S.addSpecies(5, 0, 0.5, false, &M, {"ha"});
	S.addSpecies(20, 0, 1, false, &M);
	S.resetState();
	S.initialize();
	S.print();

	return 0;
}

