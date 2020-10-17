#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model_2_ms.h"

int main(){

	TestModel M;
	Environment E;

	Solver<TestModel,Environment> S(SOLVER_FMU);
	S.addSpecies(25, 0, 1, false, &M, {}, 2);
	S.resetState();
	S.initialize();
	S.setEnvironment(&E);
	//S.print();
	
	ofstream fout("fmu_testmodel.txt");

	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out()[0] << "\t";
		//cout << S.current_time << " " [><< S.u0_out()<] << "\n";
		for (auto y : S.state) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	cout << S.u0_out()[0] << endl; 
	if (abs(S.u0_out()[0] - 1.468232) < 1e-5) return 0;
	else return 1;

}

