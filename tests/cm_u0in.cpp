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

	Solver<TestModel,Environment> S(SOLVER_CM);
	S.use_log_densities = true;
	S.addSpecies(25, 0, 1, false, &M, {}, 2);
	S.resetState();
	S.initialize();
	S.setEnvironment(&E);
	S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	//S.print();
	//for (auto s : S.state) cout << s << " "; cout << endl;

	ofstream fout("cm_testmodel.txt");

	fout << S.current_time << "\t" << 0 << "\t";
	for (auto y : S.state) fout << y << "\t"; fout << "\n";
	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out()[0] << "\t";
		//cout << S.current_time << "\t" << S.get_species(0)->size() << " " << S.u0_out()[0] << "\n";
		//cout << S.u0_out() << "\n";
		for (auto y : S.state) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	cout << S.u0_out()[0] << endl;
	// test value is from R code	
	if (abs(S.u0_out()[0] - 1.556967) < 1e-3) return 0;
	else return 1;
  
}

