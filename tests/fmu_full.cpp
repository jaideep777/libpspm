#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model.h"

int main(){

	Solver<TestModel> S(25, 0,1, SOLVER_FMU);
	
	TestModel M;
	
	S.setModel(&M);
	S.initialize();
	
	ofstream fout("fmu_testmodel.txt");

	for (double t=0; t <= 8; t=t+8.0/20) {
		S.step_to(t);
		fout << S.current_time << "\t";
		for (auto y : S.state) fout << y << "\t";
		fout << endl;
	}

	fout.close();
  
	return 0;
}

