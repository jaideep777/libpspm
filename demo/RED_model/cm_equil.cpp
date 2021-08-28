#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "red.h"

int main(){

	Species<RED_Plant> spp;
	LightEnvironment E;

	Solver S(SOLVER_CM);
	S.addSpecies(30, 1, 1e6, true, &spp, 0, 1);
	spp.set_bfin_is_u0in(false);	// say that input_birth_flux is u0
	S.use_log_densities=true;
	S.control.max_cohorts = 50;
	S.resetState();
	S.initialize();
	S.setEnvironment(&E);
	//S.print();
	
	
	ofstream fout("cm_Redmodel.txt");

	for (double t=0; t <= 5000; t=t+1) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.newborns_out()[0] << "\t";
		//cout << S.current_time << " " [><< S.u0_out()<] << "\n";
		for (auto y : S.state) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	// expected 38.1128953 (numerical R), 37.5845 (analytical)
	cout << S.newborns_out()[0] << endl; 
	//if (abs(S.u0_out()[0] - 1.468232) < 1e-5) return 0;
	//else return 1;

}

