#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model.h"

int main(){

	Solver<TestModel> S(100, 0,1, SOLVER_CM);
	
	TestModel M;
	
	S.setModel(&M);
	S.initialize();

	//for (auto s : S.state) cout << s << " "; cout << endl;

	ofstream fout("cm_testmodel.txt");

	fout << S.current_time << "\t";
	for (auto y : S.state) fout << y << "\t";
	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.newborns_out() << "\t";
		for (auto y : S.state) fout << y << "\t";
		fout << endl;
	}

	fout.close();
  
	return 0;
}

