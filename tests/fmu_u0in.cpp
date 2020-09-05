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
	S.setInputNewbornDensity(2);
	
	ofstream fout("fmu_testmodel.txt");

	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.newborns_out() << "\t";
		//cout << S.current_time << " " << S.u0_out() << "\n";
		for (auto y : S.state) fout << y << "\t";
		fout << endl;
	}

	fout.close();
 
	if (abs(S.u0_out() - 1.46823) < 1e-5) return 0;
	else return 1;

}

