#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "solver.h"
using namespace std;

#include "red.h"

int main(){

	REDModel M;
	LightEnvironment E;

	Solver<REDModel,LightEnvironment> S(SOLVER_EBT);
	S.addSpecies(30, 1, 1e6, true, &M, {});
	//S.get_species(0)->set_bfin_is_u0in(true);	// say that input_birth_flux is u0
	S.resetState();
	S.initialize();
	S.setEnvironment(&E);
	//S.print();
	
	ofstream fout("ebt_Redmodel.txt");

	for (double t=0.05; t <= 5000; t=t+100) {
		S.step_to(t);
		fout << S.current_time << "\t" << S.u0_out()[0] << "\t";
		//cout << S.current_time << " " [><< S.u0_out()<] << "\n";
		for (auto y : S.state) fout << y << "\t";
		fout << endl;
	}

	fout.close();

	cout << S.u0_out()[0] << endl; 
	//if (abs(S.u0_out()[0] - 1.468232) < 1e-5) return 0;
	//else return 1;

}

