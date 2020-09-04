#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model.h"

int main(){

	Solver<TestModel> S(25, 0,1, SOLVER_EBT);
	vector<double> x = S.getx();
	
	TestModel M;
	
	S.setModel(&M);
	S.initialize();
	
	ofstream fout("ebt_testmodel.txt");

	for (double t=0; t <= 12; t=t+0.05) {
		S.step_to(t);
	
		vector<double> v = S.cohortsToDensity_EBT(x);
		
		fout << S.current_time << "\t";
		for (auto y : S.state) fout << y << "\t";
		fout << endl;
	}

	fout.close();
  
	return 0;
}

