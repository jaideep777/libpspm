#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "test_model_2_ms.h"

std::vector <double> myseq(double from, double to, int len){
	std::vector<double> x(len);
	for (size_t i=0; i<len; ++i) x[i] = from + i*(to-from)/(len-1);
	return x;
}

int main(){

	Species<TestModel> spp;
	Environment E;

	Solver S(SOLVER_CM);
	S.use_log_densities = true;
	S.control.cm_grad_dx = 0.001;
	S.control.max_cohorts = 26;
	S.addSpecies(25, 0, 1, false, &spp, 4, -1);
	S.resetState();
	S.initialize();
	S.setEnvironment(&E);
	S.species_vec[0]->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	S.print();
	//for (auto s : S.state) cout << s << " "; cout << endl;

	ofstream fout("cm_testmodel_equil.txt");

	//fout << S.current_time << "\t" << 0 << "\t";
	//for (auto y : S.state){fout << y << "\t";} fout << "\n";
	for (double t=0.05; t <= 8; t=t+0.05) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out(t)[0] << "\t";
		cout << S.current_time << "\t" << S.u0_out(t)[0] << "\n";
		//cout << S.current_time << "\t" << S.species_vec[0]->xsize() << " " << S.u0_out()[0] << "\t" << S.species_vec[0]->get_boundary_u() << "\n";
		//cout << S.u0_out() << "\n";
		vector<double> breaks = myseq(0,1,26);
		vector<double> v = S.getDensitySpecies(0, breaks);
		for (auto y : v) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	cout << S.u0_out(S.current_time)[0] << endl;
	// test value is from R code	
	//if (abs(S.u0_out()[0] - 1.556967) < 1e-5) return 0;  // this is when integrate_x BC is not included
	if (abs(S.u0_out(S.current_time)[0] - 0.976177) < 1e-5) return 0;  // this is when integrate_x BC IS included

	else return 1;
  
}
